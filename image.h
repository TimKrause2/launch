#ifndef IMAGE_H
#define IMAGE_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	unsigned char r;
	unsigned char g;
	unsigned char b;
} pixel;

typedef struct {
    unsigned char r;
    unsigned char g;
    unsigned char b;
    unsigned char a;
} rgba;

typedef struct {
	int width;
	int height;
	pixel *data;
} image;

image * image_load( const char *filename );
image * image_new( int width, int height );

void image_free( image *arg );

#ifdef __cplusplus
}
#endif

#endif
