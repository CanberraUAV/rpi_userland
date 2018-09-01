/*
  extract RPI raw10 image, producing a 16 bit pgm and 8 bit ppm
  
  With thanks to http://github.com/6by9/RPiTest
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>
#include <sys/wait.h>
#include <signal.h>       
#include <math.h>       

#include "libjpeg/jpeglib.h"
#include <stdbool.h>
#include "cuav_util.h"
#include <sys/time.h>
#include <sys/mman.h>
#include <pthread.h>

#pragma GCC optimize("O3")

// offset from end of the to BRCM marker
#define BRCM_OFFSET 10270208
#define DATA_OFFSET 0x8000

// RPI image size
#define SCALING 2
#define IMG_WIDTH (3280/SCALING)
#define IMG_HEIGHT (2464/SCALING)

static unsigned num_children_created;
static pthread_mutex_t counter_lock = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t jpeg_lock = PTHREAD_MUTEX_INITIALIZER;
static volatile unsigned num_children_exited;

static void thread_exit(void)
{
    pthread_mutex_lock(&counter_lock);
    num_children_exited++;
    pthread_mutex_unlock(&counter_lock);
    static int ret = 0;
    pthread_exit(&ret);
}
    

#define PACKED __attribute__((__packed__))

struct PACKED rgb8 {
    uint8_t r, g, b;
};

struct PACKED rgbf {
    float r, g, b;
};

/*
  16 bit bayer grid
 */
struct PACKED bayer_image {
    uint16_t data[IMG_HEIGHT*SCALING][IMG_WIDTH*SCALING];
};

/*
  RGB image, 8 bit
 */
struct PACKED rgb8_image {
    struct rgb8 data[IMG_HEIGHT][IMG_WIDTH];
};

/*
  RGB image, float
 */
struct PACKED rgbf_image {
    struct rgbf data[IMG_HEIGHT][IMG_WIDTH];
};

/*
  len shading scaling array
 */
struct lens_shading {
    float scale[IMG_HEIGHT][IMG_WIDTH];
};

struct brcm_header {
    char tag[4]; // BRCM
    uint8_t pad[172];
    uint8_t name[32];
    uint16_t width;
    uint16_t height;
    uint16_t padding_right;
    uint16_t padding_down;
    uint32_t dummy[6];
    uint16_t transform;
    uint16_t format;
    uint8_t bayer_order;
    uint8_t bayer_format;
};

extern void swab(const void *from, void *to, ssize_t n);

static void extract_raw10(const uint8_t *b, uint16_t width, uint16_t height, uint16_t raw_stride, struct bayer_image *bayer)
{
    uint8_t data[raw_stride];
    uint16_t row, col;
    
    for (row=0; row<height; row++) {
        uint16_t *raw = &bayer->data[row][0];
        memcpy(data, b, raw_stride);
        b += raw_stride;
        uint8_t *dp = &data[0];
        for (col=0; col<width; col+=4, dp+=5) {
            // the top two bits are packed into byte 4 of each group
            raw[col+0] = dp[0] << 2 | (dp[4]&3);
            raw[col+1] = dp[1] << 2 | ((dp[4]>>2)&3);
            raw[col+2] = dp[2] << 2 | ((dp[4]>>4)&3);
            raw[col+3] = dp[3] << 2 | ((dp[4]>>6)&3);
        }
    }
}


#if SCALING == 1
static void debayer_BGGR_float(const struct bayer_image *bayer, struct rgbf_image *rgb)
{
    /*
      layout in the input image is in blocks of 4 values. The top
      left corner of the image looks like this
      B G B G
      G R G R
      B G B G
      G R G R
    */
    uint16_t x, y;
    for (y=1; y<IMG_HEIGHT-2; y += 2) {
        for (x=1; x<IMG_WIDTH-2; x += 2) {
            rgb->data[y+0][x+0].r = bayer->data[y+0][x+0];
            rgb->data[y+0][x+0].g = ((uint16_t)bayer->data[y-1][x+0] + (uint16_t)bayer->data[y+0][x-1] +
                                     (uint16_t)bayer->data[y+1][x+0] + (uint16_t)bayer->data[y+0][x+1]) >> 2;
            rgb->data[y+0][x+0].b = ((uint16_t)bayer->data[y-1][x-1] + (uint16_t)bayer->data[y+1][x-1] +
                                     (uint16_t)bayer->data[y-1][x+1] + (uint16_t)bayer->data[y+1][x+1]) >> 2;
            rgb->data[y+0][x+1].r = ((uint16_t)bayer->data[y+0][x+0] + (uint16_t)bayer->data[y+0][x+2]) >> 1;
            rgb->data[y+0][x+1].g = bayer->data[y+0][x+1];
            rgb->data[y+0][x+1].b = ((uint16_t)bayer->data[y-1][x+1] + (uint16_t)bayer->data[y+1][x+1]) >> 1;

            rgb->data[y+1][x+0].r = ((uint16_t)bayer->data[y+0][x+0] + (uint16_t)bayer->data[y+2][x+0]) >> 1;
            rgb->data[y+1][x+0].g = bayer->data[y+1][x+0];
            rgb->data[y+1][x+0].b = ((uint16_t)bayer->data[y+1][x-1] + (uint16_t)bayer->data[y+1][x+1]) >> 1;

            rgb->data[y+1][x+1].r = ((uint16_t)bayer->data[y+0][x+0] + (uint16_t)bayer->data[y+2][x+0] +
                                     (uint16_t)bayer->data[y+0][x+2] + (uint16_t)bayer->data[y+2][x+2]) >> 2;
            rgb->data[y+1][x+1].g = ((uint16_t)bayer->data[y+0][x+1] + (uint16_t)bayer->data[y+1][x+2] +
                                     (uint16_t)bayer->data[y+2][x+1] + (uint16_t)bayer->data[y+1][x+0]) >> 2;
            rgb->data[y+1][x+1].b = bayer->data[y+1][x+1];
        }
        rgb->data[y+0][0] = rgb->data[y+0][1];
        rgb->data[y+1][0] = rgb->data[y+1][1];
        rgb->data[y+0][IMG_WIDTH-1] = rgb->data[y+0][IMG_WIDTH-2];
        rgb->data[y+1][IMG_WIDTH-1] = rgb->data[y+1][IMG_WIDTH-2];
    }
    memcpy(rgb->data[0], rgb->data[1], IMG_WIDTH*sizeof(rgb->data[0][0]));
    memcpy(rgb->data[IMG_HEIGHT-1], rgb->data[IMG_HEIGHT-2], IMG_WIDTH*sizeof(rgb->data[0][0]));
}
#elif SCALING == 2
static void debayer_BGGR_float(const struct bayer_image *bayer, struct rgbf_image *rgb)
{
    /*
      layout in the input image is in blocks of 4 values. The top
      left corner of the image looks like this
      B G B G
      G R G R
      B G B G
      G R G R
    */
    uint16_t x, y;
    for (y=1; y<IMG_HEIGHT*SCALING-2; y += 2) {
        for (x=1; x<IMG_WIDTH*SCALING-2; x += 2) {
            float r, g, b;
            
            r = bayer->data[y+0][x+0];
            g = ((uint16_t)bayer->data[y-1][x+0] + (uint16_t)bayer->data[y+0][x-1] +
                 (uint16_t)bayer->data[y+1][x+0] + (uint16_t)bayer->data[y+0][x+1]) >> 2;
            b = ((uint16_t)bayer->data[y-1][x-1] + (uint16_t)bayer->data[y+1][x-1] +
                 (uint16_t)bayer->data[y-1][x+1] + (uint16_t)bayer->data[y+1][x+1]) >> 2;
            
            r += ((uint16_t)bayer->data[y+0][x+0] + (uint16_t)bayer->data[y+0][x+2]) >> 1;
            g += bayer->data[y+0][x+1];
            b += ((uint16_t)bayer->data[y-1][x+1] + (uint16_t)bayer->data[y+1][x+1]) >> 1;

            r += ((uint16_t)bayer->data[y+0][x+0] + (uint16_t)bayer->data[y+2][x+0]) >> 1;
            g += bayer->data[y+1][x+0];
            b += ((uint16_t)bayer->data[y+1][x-1] + (uint16_t)bayer->data[y+1][x+1]) >> 1;

            r += ((uint16_t)bayer->data[y+0][x+0] + (uint16_t)bayer->data[y+2][x+0] +
                  (uint16_t)bayer->data[y+0][x+2] + (uint16_t)bayer->data[y+2][x+2]) >> 2;
            g += ((uint16_t)bayer->data[y+0][x+1] + (uint16_t)bayer->data[y+1][x+2] +
                  (uint16_t)bayer->data[y+2][x+1] + (uint16_t)bayer->data[y+1][x+0]) >> 2;
            b += bayer->data[y+1][x+1];
            
            rgb->data[y/2][x/2].r = r*0.25;
            rgb->data[y/2][x/2].g = g*0.25;
            rgb->data[y/2][x/2].b = b*0.25;
        }
        //rgb->data[y/2][0] = rgb->data[y/2][1];
        //rgb->data[y+1][0] = rgb->data[y+1][1];
        //rgb->data[y+0][IMG_WIDTH-1] = rgb->data[y+0][IMG_WIDTH-2];
        //rgb->data[y+1][IMG_WIDTH-1] = rgb->data[y+1][IMG_WIDTH-2];
    }
    //memcpy(rgb->data[0], rgb->data[1], IMG_WIDTH*sizeof(rgb->data[0][0]));
    //memcpy(rgb->data[IMG_HEIGHT-1], rgb->data[IMG_HEIGHT-2], IMG_WIDTH*sizeof(rgb->data[0][0]));
}
#endif

static struct lens_shading *shading;
const float shading_scale_factor = 1.9;

/*
  create a lens shading correction array
 */
static void create_lens_shading(void)
{
    shading = malloc(sizeof(*shading));
    uint16_t y, x;
    for (y=0; y<IMG_HEIGHT; y++) {
        for (x=0; x<IMG_WIDTH; x++) {
            float dx = fabsf(((float)x) - IMG_WIDTH/2) / (IMG_WIDTH/2);
            float dy = fabsf(((float)y) - IMG_HEIGHT/2) / (IMG_HEIGHT/2);
            float from_center = sqrt(dx*dx + dy*dy);
            if (from_center > 1.0) {
                from_center = 1.0;
            }
            shading->scale[y][x] = 1.0 + from_center * shading_scale_factor;
        }
    }
}

static void rgbf_to_rgb8(const struct rgbf_image *rgbf, struct rgb8_image *rgb8)
{
    float highest = 0;
    uint16_t x, y;
    
    for (y=0; y<IMG_HEIGHT; y++) {
        for (x=0; x<IMG_WIDTH; x++) {
            const struct rgbf *d = &rgbf->data[y][x];
            if (d->r > highest) {
                highest = d->r * shading->scale[y][x];
            }
            if (d->g > highest) {
                highest = d->g * shading->scale[y][x];
            }
            if (d->b > highest) {
                highest = d->b * shading->scale[y][x];
            }
        }
    }

    float scale = 255 / highest;
    const float cscale[3] = { 1, 0.48, 0.82 };
#define MIN(a,b) ((a)<(b)?(a):(b))
    for (y=0; y<IMG_HEIGHT; y++) {
        for (x=0; x<IMG_WIDTH; x++) {
            float shade_scale = shading->scale[y][x];
            if (rgbf->data[y][x].r >= 1022) {
                rgb8->data[y][x].r = 255;
            } else {
                rgb8->data[y][x].r = MIN(rgbf->data[y][x].r * scale * cscale[0] * shade_scale, 255);
            }
            if (rgbf->data[y][x].g >= 1022) {
                rgb8->data[y][x].g = 255;
            } else {
                rgb8->data[y][x].g = MIN(rgbf->data[y][x].g * scale * cscale[1] * shade_scale, 255);
            }
            if (rgbf->data[y][x].b >= 1022) {
                rgb8->data[y][x].b = 255;
            } else {
                rgb8->data[y][x].b = MIN(rgbf->data[y][x].b * scale * cscale[2] * shade_scale, 255);
            }
        }
    }
}

/*
  extract bayer data from a RPi image
 */
static void extract_rpi_bayer(const uint8_t *buffer, uint32_t size, struct bayer_image *bayer)
{
    const uint8_t *b;
    b = &buffer[size-BRCM_OFFSET];
    struct brcm_header header;

    memcpy(&header, b, sizeof(header));

    if (strncmp(header.tag, "BRCM", 4) != 0) {
        printf("bad header name - expected BRCM\n");
        thread_exit();
    }

    uint32_t raw_stride = ((((((header.width + header.padding_right)*5)+3)>>2) + 31)&(~31));
    
    printf("Image %ux%u format %u '%s' stride:%u bayer_order:%u\n",
           header.width, header.height, header.format, header.name, raw_stride,
           header.bayer_order);

    if (header.width != IMG_WIDTH*SCALING || header.height != IMG_HEIGHT*SCALING) {
        printf("Unexpected image size\n");
        thread_exit();
    }

    b = &buffer[size - (BRCM_OFFSET - DATA_OFFSET)];
    
    extract_raw10(b, header.width, header.height, raw_stride, bayer);
}

/*
  write a JPG image
 */
static bool write_JPG(const char *filename, const struct rgb8_image *img, int quality)
{
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE *outfile;
    
    pthread_mutex_lock(&jpeg_lock);
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);

    if ((outfile = fopen(filename, "wb")) == NULL) {
        fprintf(stderr, "can't open %s\n", filename);
        pthread_mutex_unlock(&jpeg_lock);
        return false;
    }
    jpeg_stdio_dest(&cinfo, outfile);

    cinfo.image_width = IMG_WIDTH;
    cinfo.image_height = IMG_HEIGHT;
    cinfo.input_components = 3;
    cinfo.in_color_space = JCS_RGB;

    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE);

    jpeg_start_compress(&cinfo, TRUE);

    while (cinfo.next_scanline < cinfo.image_height) {
        JSAMPROW row[1];
        row[0] = (JSAMPROW)&img->data[cinfo.next_scanline][0];
        jpeg_write_scanlines(&cinfo, row, 1);
    }
    
    jpeg_finish_compress(&cinfo);
    pthread_mutex_unlock(&jpeg_lock);
    fclose(outfile);
    jpeg_destroy_compress(&cinfo);

    return true;    
}

/*
  automatically cope with system load
 */
static void control_delay(void)
{
    static unsigned delay_us = 100000;
    pthread_mutex_lock(&counter_lock);
    int children_active = (int)num_children_created - (int)num_children_exited;
    pthread_mutex_unlock(&counter_lock);
    if (children_active > 12) {
        delay_us *= 1.2;
    } else if (children_active < 8) {
        delay_us *= 0.9;
    }
    if (delay_us < 1000) {
        delay_us = 1000;
    }
    printf("Delay %u active %d\n", delay_us, children_active);
    usleep(delay_us);
}

struct child_info {
    char *fname;
    const char *linkname;
    uint8_t *buffer;
    uint32_t size;
};

static void *process_image(void *arg)
{
    struct child_info *cinfo = arg;

    struct bayer_image *bayer = mm_alloc(sizeof(*bayer));
    if (!bayer) {
        mm_free(cinfo->buffer, cinfo->size);
        free(cinfo->fname);
        free(cinfo);
        thread_exit();
        return NULL;
    }
    
    extract_rpi_bayer(cinfo->buffer, cinfo->size, bayer);

    mm_free((void *)cinfo->buffer, cinfo->size);
        
    struct rgbf_image *rgbf = mm_alloc(sizeof(*rgbf));
    if (!rgbf) {
        mm_free(bayer, sizeof(*bayer));
        free(cinfo->fname);
        free(cinfo);
        thread_exit();
        return NULL;
    }
    
    debayer_BGGR_float(bayer, rgbf);

    mm_free(bayer, sizeof(*bayer));

    struct rgb8_image *rgb8 = mm_alloc(sizeof(*rgb8));
    if (!rgb8) {
        mm_free(rgbf, sizeof(*rgbf));
        free(cinfo->fname);
        free(cinfo);
        thread_exit();
        return NULL;
    }
    
    rgbf_to_rgb8(rgbf, rgb8);

    mm_free(rgbf, sizeof(*rgbf));
    
    char *stop_name = NULL;
    struct stat st;
    asprintf(&stop_name, "%s.stop", cinfo->linkname);
    if (stop_name && stat(stop_name, &st) == 0) {
        printf("** Stop file exists\n");
        thread_exit();
        return NULL;
    }
    write_JPG(cinfo->fname, rgb8, 100);
    unlink(cinfo->linkname);
    symlink(cinfo->fname, cinfo->linkname);

    mm_free(rgb8, sizeof(*rgb8));
    free(cinfo->fname);
    free(cinfo);
    thread_exit();
    return NULL;
}

void cuav_process(const uint8_t *buffer, uint32_t size, const char *filename, const char *linkname, const struct timeval *tv, bool halfres)
{
    printf("Processing %u bytes\n", size);

    struct tm tm;
    time_t t = tv->tv_sec;
    gmtime_r(&t, &tm);

    if (!shading) {
        create_lens_shading();
    }

    char *fname = NULL;
    asprintf(&fname, "%s%04u%02u%02u%02u%02u%02u%02uZ.jpg",
             filename,
             tm.tm_year+1900,
             tm.tm_mon+1,
             tm.tm_mday,
             tm.tm_hour,
             tm.tm_min,
             tm.tm_sec,
             (unsigned)(tv->tv_usec/10000));
    if (!fname) {
        return;
    }
    printf("fname=%s\n", fname);

    num_children_created++;

    uint8_t *buf2 = mm_alloc(size);
    if (!buf2) {
        free(fname);
        return;
    }
    memcpy(buf2, buffer, size);

    struct child_info *cinfo;
    cinfo = malloc(sizeof(*cinfo));
    if (!cinfo) {
        mm_free(buf2, size);
        free(fname);
        return;
    }
    cinfo->fname = fname;
    cinfo->buffer = buf2;
    cinfo->size = size;
    cinfo->linkname = linkname;

    pthread_t thread;
    if (pthread_create(&thread, NULL, process_image, cinfo) != 0) {
        mm_free(buf2, size);
        free(fname);
        free(cinfo);
        return;        
    }
    pthread_detach(thread);

    control_delay();
}

void *mm_alloc(uint32_t size)
{
    uint32_t pagesize = getpagesize();
    uint32_t num_pages = (size + pagesize - 1) / pagesize;
    return mmap(0, num_pages*pagesize, PROT_READ | PROT_WRITE, 
                MAP_ANON | MAP_PRIVATE, -1, 0);
}

void mm_free(void *ptr, uint32_t size)
{
    uint32_t pagesize = getpagesize();
    uint32_t num_pages = (size + pagesize - 1) / pagesize;
    munmap(ptr, num_pages*pagesize);
}
