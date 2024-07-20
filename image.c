#include <stdio.h>
#include <jpeglib.h>
#include <tiffio.h>
#include "image.h"

#include <stdlib.h>
#include <string.h>

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

void image_free( image *arg )
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




