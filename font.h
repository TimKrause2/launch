//#include <GL/gl.h>
#include <GLES/gl32.h>
#include <GLES3/gl3ext.h>

#include <vector>

#ifdef _MSC_VER
#define MIN __min
#define MAX __max
#else
#define MIN std::min
#define MAX std::max
#endif


// Define some fixed size types.

typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;


// A simple 32-bit pixel.

union Pixel32
{
  Pixel32()
  : integer(0) { }
  Pixel32(uint8 bi, uint8 gi, uint8 ri, uint8 ai = 255)
  {
    b = bi;
    g = gi;
    r = ri;
    a = ai;
  }

  uint32 integer;

  struct
  {
    uint8 r, g, b, a;
  };
};


struct Vec2
{
  Vec2() { }
  Vec2(float a, float b)
  : x(a), y(b) { }

  float x, y;
};


struct Rect
{
  Rect() { }
  Rect(float left, float top, float right, float bottom)
  : xmin(left), xmax(right), ymin(top), ymax(bottom) { }

  void Include(const Vec2 &r)
  {
    xmin = MIN(xmin, r.x);
    ymin = MIN(ymin, r.y);
    xmax = MAX(xmax, r.x);
    ymax = MAX(ymax, r.y);
  }

  float Width() const { return xmax - xmin + 1; }
  float Height() const { return ymax - ymin + 1; }

  float xmin, xmax, ymin, ymax;
};

// A horizontal pixel span generated by the FreeType renderer.

struct Span
{
  Span() { }
  Span(int _x, int _y, int _width, int _coverage)
  : x(_x), y(_y), width(_width), coverage(_coverage) { }

  int x, y, width, coverage;
};

typedef std::vector<Span> Spans;

class Font
{
private:
	GLuint *textures;
	GLuint list_base;
	void Free( void );
public:
	Font( void );
	~Font( void );
	
	void Load( const char *filename, unsigned int pixelsize );
	void LoadOutline( const char *filename, unsigned int fontsize,
										Pixel32 fontColor, Pixel32 outlineColor,
									 double outlineWidth);
	void Printf( double x, double y, const char *format, ... );
};