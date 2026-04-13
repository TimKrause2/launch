#include <stdio.h>
#include <jpeglib.h>
#include <tiffio.h>
#include "image.h"
#include "shader.hpp"
#include <stdlib.h>
#include <string.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
//#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>

static GLuint program;
static GLuint vao; // vertex array object
static GLuint vbuff;
static GLuint uvbuff;
static GLuint mvp_id;
static GLuint tex_id;
static GLuint alpha_id;

void image_init(void)
{
    program = LoadShaders( "TransformVertexShader.vertexshader",
                           "image.fgmt.glsl" );
    mvp_id = glGetUniformLocation(program, "MVP");
    tex_id  = glGetUniformLocation(program, "myTextureSampler");
    alpha_id = glGetUniformLocation(program, "alpha");

    glm::vec3 vertices[4] =
    {
        {-0.5f, 0.5f, 0.0f},
        {-0.5f,-0.5f, 0.0f},
        { 0.5f,-0.5f, 0.0f},
        { 0.5f, 0.5f, 0.0f}
    };

    glm::vec2 uvs[4] =
    {
        {0.0f, 0.0f},
        {0.0f, 1.0f},
        {1.0f, 1.0f},
        {1.0f, 0.0f}
    };

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glGenBuffers(1, &vbuff);
    glBindBuffer(GL_ARRAY_BUFFER, vbuff);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glGenBuffers(1, &uvbuff);
    glBindBuffer(GL_ARRAY_BUFFER, uvbuff);
    glBufferData(GL_ARRAY_BUFFER, sizeof(uvs), uvs, GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vbuff);
    glVertexAttribPointer(
        0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
        3,                  // size
        GL_FLOAT,           // type
        GL_FALSE,           // normalized?
        0,                  // stride
        (void*)0            // array buffer offset
    );

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, uvbuff);
    glVertexAttribPointer(
        1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
        2,                                // size : U+V => 2
        GL_FLOAT,                         // type
        GL_FALSE,                         // normalized?
        0,                                // stride
        (void*)0                          // array buffer offset
    );

    glBindVertexArray(0);
}

glTexture* texture_load( const char *tex_file )
{
    glTexture *result = (glTexture*)malloc(sizeof(glTexture));

    // load the texture for the body
    image* tex_image = image_load(tex_file);
    if(!tex_image){
        printf("couldn't load texture file:\"%s\"",tex_file);
        return 0;
    }

    GLuint texName;

    glGenTextures(1, &texName);
    glBindTexture(GL_TEXTURE_2D, texName);
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );

    glTexImage2D(
                GL_TEXTURE_2D, 0, GL_RGB8,
                tex_image->width, tex_image->height,
                0,
                GL_RGB, GL_UNSIGNED_BYTE,
                (const GLvoid*)tex_image->data );

    result->width = tex_image->width;
    result->height = tex_image->height;
    result->tex_id = texName;

    image_free(tex_image);

    return result;
}

glTexture* rgba_texture_load( const char *tex_file )
{
    GLuint texName;
    glTexture *result = (glTexture*)malloc(sizeof(glTexture));

    // load the texture for the body
    rgba_image* tex_image = rgba_image_load(tex_file);
    if(!tex_image){
        printf("couldn't load texture file:\"%s\"",tex_file);
        return 0;
    }

    glGenTextures(1, &texName);
    glBindTexture(GL_TEXTURE_2D, texName);
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );

    glTexImage2D(
                GL_TEXTURE_2D, 0, GL_RGBA8,
                tex_image->width, tex_image->height,
                0,
                GL_RGBA, GL_UNSIGNED_BYTE,
                (const GLvoid*)tex_image->data );

    result->width = tex_image->width;
    result->height = tex_image->height;
    result->tex_id = texName;

    rgba_image_free(tex_image);


    return result;
}

void texture_free( glTexture* tex )
{
    glDeleteTextures(1, &tex->tex_id);
    free(tex);
}

void texture_sprite(glTexture *tex, float xc, float yc, float alpha)
{
    // initialize the scale transform
    glm::vec3 vscale(tex->width, tex->height, 1.0f);
    glm::mat4 Mscale = glm::scale(vscale);

    // initialize the translation matrix
    glm::vec3 vtrans(xc, yc, 0.0f);
    glm::mat4 Mtrans = glm::translate(vtrans);

    // initialize the projection matrix
    GLint vp[4];
    glGetIntegerv(GL_VIEWPORT, vp);
    float left = 0.0f;
    float right = vp[2];
    float bottom = 0.0f;
    float top = vp[3];
    float near = 1.0f;
    float far = -1.0f;
    glm::mat4 Mproj = glm::ortho(left, right, bottom, top,
                                 near, far);

    // calculate the mvp matrix
    glm::mat4 mvp = Mproj*Mtrans*Mscale;

    glBindVertexArray(vao);
    glUseProgram(program);
    glUniformMatrix4fv(mvp_id, 1, GL_FALSE, glm::value_ptr(mvp));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex->tex_id);
    glUniform1i(tex_id, 0);

    glUniform1f(alpha_id, alpha);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);

    glDisable(GL_BLEND);
    glBindVertexArray(0);
    glUseProgram(0);
}

image * image_new( int width, int height )
{
    image *r = (image*)malloc(sizeof(image));
    if(!r) return NULL;
    r->data = (pixel*)malloc(sizeof(pixel)*width*height);
    if(!r->data){
        free(r);
        return NULL;
    }
    r->width = width;
    r->height = height;
    return r;
}

rgba_image * rgba_image_new( int width, int height )
{
    rgba_image *r = (rgba_image*)malloc(sizeof(rgba_image));
    if(!r) return NULL;
    r->data = (rgba*)malloc(sizeof(rgba)*width*height);
    if(!r->data){
        free(r);
        return NULL;
    }
    r->width = width;
    r->height = height;
    return r;
}

void image_free( image *arg )
{
    free( arg->data );
    free( arg );
}

void rgba_image_free( rgba_image *arg )
{
    free( arg->data );
    free( arg );
}

const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

image* image_load_jpeg(const char *filename);
image* image_load_tiff(const char *filename);
rgba_image* rgba_image_load_jpeg(const char *filename);
rgba_image* rgba_image_load_tiff(const char *filename);

image* image_load(const char *filename)
{
    const char *ext = get_filename_ext(filename);
    if(!*ext){
        printf("Can't determine file type.\n");
        return NULL;
    }
    if(strcmp(ext,"jpg")==0){
        return image_load_jpeg(filename);
    }
    if(strcmp(ext,"tif")==0 || strcmp(ext,"tiff")==0){
        return image_load_tiff(filename);
    }
    return NULL;
}

rgba_image* rgba_image_load(const char *filename)
{
    const char *ext = get_filename_ext(filename);
    if(!*ext){
        printf("Can't determine file type.\n");
        return NULL;
    }
    if(strcmp(ext,"jpg")==0){
        return rgba_image_load_jpeg(filename);
    }
    if(strcmp(ext,"tif")==0 || strcmp(ext,"tiff")==0){
        return rgba_image_load_tiff(filename);
    }
    return NULL;
}

image* image_load_tiff(const char *filename){
    image* result;
    TIFF* tif = TIFFOpen(filename, "r");
    if(tif){
        uint32 w, h;
        size_t npixels;
        uint32* rgba_data;

        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
        result = image_new(w, h);
        if(!result){
            printf("Couldn't allocate image\n");
            TIFFClose(tif);
            return NULL;
        }
        npixels = w * h;
        rgba_data = (uint32*) _TIFFmalloc(npixels*sizeof(uint32));
        if(rgba_data){
            if(TIFFReadRGBAImage(tif, w, h, rgba_data, 0)){
                // convert rgba to the image rgb
                rgba* src = (rgba*)rgba_data;
                pixel* dst = result->data;
                for(int y=0;y<h;y++){
                    for(int x=0;x<w;x++){
                        dst->r = src->r;
                        dst->g = src->g;
                        dst->b = src->b;
                        src++;
                        dst++;
                    }
                }
            }
            _TIFFfree(rgba_data);
        }else{
            printf("Couldn't allocate tiff data\n");
            image_free(result);
            TIFFClose(tif);
            return NULL;
        }

        TIFFClose(tif);

        return result;
    }else{
        printf("Couldn't open tiff \"%s\"\n",filename);
        return NULL;
    }
    return NULL;
}

rgba_image* rgba_image_load_tiff(const char *filename){
    rgba_image* result;
    TIFF* tif = TIFFOpen(filename, "r");
    if(tif){
        uint32 w, h;
        size_t npixels;
        uint32* rgba_data;

        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
        result = rgba_image_new(w, h);
        if(!result){
            printf("Couldn't allocate image\n");
            TIFFClose(tif);
            return NULL;
        }
        npixels = w * h;
        rgba_data = (uint32*) _TIFFmalloc(npixels*sizeof(uint32));
        if(rgba_data){
            if(TIFFReadRGBAImage(tif, w, h, rgba_data, 0)){
                // convert rgba to the image rgb
                rgba* src = (rgba*)rgba_data;
                rgba* dst = result->data;
                for(int y=0;y<h;y++){
                    for(int x=0;x<w;x++){
                        *dst = *src;
                        src++;
                        dst++;
                    }
                }
            }
            _TIFFfree(rgba_data);
        }else{
            printf("Couldn't allocate tiff data\n");
            rgba_image_free(result);
            TIFFClose(tif);
            return NULL;
        }

        TIFFClose(tif);

        return result;
    }else{
        printf("Couldn't open tiff \"%s\"\n",filename);
        return NULL;
    }
    return NULL;
}

image * image_load_jpeg( const char *filename )
{
    image * r;
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    FILE * infile;

    infile = fopen( filename, "rb" );
    if(!infile){
        perror("fopen");
        return NULL;
    }

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);

    jpeg_stdio_src( &cinfo, infile );

    jpeg_read_header( &cinfo, TRUE );

    cinfo.out_color_space = JCS_RGB;

    jpeg_start_decompress( &cinfo );

    printf("cinfo.out_color_space:%d\n",cinfo.out_color_space);
    printf("cinfo.output_components:%d\n",cinfo.output_components);


    r = image_new( cinfo.output_width, cinfo.output_height );
    if( !r ){
        printf("Couldn't allocate image data\n");
        jpeg_destroy_decompress( &cinfo );
        return NULL;
    }else{
        printf("Allocated image memory for \"%s\" at %dx%d pixels\n",
               filename, cinfo.output_width, cinfo.output_height);
    }

    JSAMPROW scanlines[1];
    while (cinfo.output_scanline < cinfo.output_height) {
        scanlines[0] = (JSAMPROW) &(r->data[cinfo.output_scanline*cinfo.output_width]);
        jpeg_read_scanlines( &cinfo, scanlines, 1 );
// 		pixel *p = (pixel*)scanlines[0];
// 		for(int i=0;i<cinfo.output_width;i++,p++){
// 			printf("%02X%02X%02X",p->r,p->g,p->b);
// 		}
// 		printf("\n");
    }

    jpeg_finish_decompress( &cinfo );

    jpeg_destroy_decompress( &cinfo );

    fclose( infile );


    return r;
}

rgba_image * rgba_image_load_jpeg( const char *filename )
{
    rgba_image * r;
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    FILE * infile;

    infile = fopen( filename, "rb" );
    if(!infile){
        perror("fopen");
        return NULL;
    }

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);

    jpeg_stdio_src( &cinfo, infile );

    jpeg_read_header( &cinfo, TRUE );

    cinfo.out_color_space = JCS_RGB;

    jpeg_start_decompress( &cinfo );

    printf("cinfo.out_color_space:%d\n",cinfo.out_color_space);
    printf("cinfo.output_components:%d\n",cinfo.output_components);


    r = rgba_image_new( cinfo.output_width, cinfo.output_height );
    if( !r ){
        printf("Couldn't allocate image data\n");
        jpeg_destroy_decompress( &cinfo );
        return NULL;
    }else{
        printf("Allocated image memory for \"%s\" at %dx%d pixels\n",
               filename, cinfo.output_width, cinfo.output_height);
    }

    pixel *pixels = (pixel*)malloc(sizeof(pixel)*cinfo.output_width);
    if(!pixels){
        printf("Couldn't allocate scanline storage.\n");
    }

    JSAMPROW scanlines[1];
    rgba *dst = r->data;
    while (cinfo.output_scanline < cinfo.output_height) {
        scanlines[0] = (JSAMPROW) pixels;
        jpeg_read_scanlines( &cinfo, scanlines, 1 );
        pixel *p = pixels;
        for(int i=0;i<cinfo.output_width;i++){
            dst->r = p->r;
            dst->g = p->g;
            dst->b = p->b;
            dst->a = 255;
            dst++;
            p++;
        }
    }

    free(pixels);

    jpeg_finish_decompress( &cinfo );

    jpeg_destroy_decompress( &cinfo );

    fclose( infile );

    return r;
}




