#ifndef IMAGE_H
#define IMAGE_H

#include <GLES3/gl32.h>
#include <GLES3/gl3ext.h>

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

typedef struct {
    int width;
    int height;
    rgba *data;
} rgba_image;

typedef struct {
    int width;
    int height;
    GLuint tex_id;
} glTexture;

void image_init(void);


glTexture* texture_load( const char *filename );
glTexture* rgba_texture_load( const char *filename );
void texture_free( glTexture* tex);
void texture_sprite(glTexture *tex,
                    float xc,
                    float yc, // origin is bottom left corner. positive yc moves up
                    float alpha);

image * image_load( const char *filename );
image * image_new( int width, int height );

rgba_image* rgba_image_load( const char *filename );
rgba_image* rgba_image_new( int width, int height );

void image_free( image *arg );
void rgba_image_free( rgba_image *img );

#ifdef __cplusplus
}
#endif

#endif
