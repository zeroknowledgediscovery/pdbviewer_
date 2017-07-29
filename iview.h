#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GL/gl.h>
#include <GL/gle.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <jpeglib.h>
#include <boost/algorithm/string.hpp>
#include "config.h"
#include <memory.h>
#include <sys/time.h>
#include <gl2ps.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <exception>
#include <boost/timer/timer.hpp>

#define X_AXIS                  0
#define Y_AXIS                  1
#define Z_AXIS                  2

#define DEBUG_                  1

using namespace boost::algorithm;
using namespace std;
using namespace boost::program_options;
//-------------------------------------------------

/* vector */
typedef float vec3_t[3]; 

/* enumerations for the mouse buttons */
enum {
  UP = 1, 
  DOWN, 
};
enum {
  CAPTURE_VIEWPORT=1,
  ATOM_CA,
  ATOM_CA_CB,
  ATOM_CA_CB_N,
  ATOM_CA_CB_N_O,
  ANIMATE,
  MENU_EXIT
};


/*!   read config file and set up parameters. */
void molviewConfig(string specstr="");


/*! Write JPEG file from sceenshot  using libjpeg */
void write_JPEG_file (const char *filename, 
		      int quality,
		      int image_width, 
		      int image_height);


/*! Function to print screen */
void CaptureViewPort() ;

/*! Function to print screen in high res*/
void CaptureViewPort_HR();

/*!   glInit   Sets up various OpenGL stuff. */
void glInit();

void update_show(vector<string> str);
void SelectFromMenu(int idCommand);

int BuildPopupMenu (void);
//-------------------------------------------------
static double vlen(double,double,double);
static void   pos(double*,double*,double*,const int,const int,const int*);
static void   getMatrix();
static void   invertMatrix(const GLdouble *m, GLdouble *out );

void glutMotion(int, int) ;
void glutMouse(int, int, int, int) ;
void zprSelectionFunc(void (*f)(void));
void zprPickFunc(void (*f)(GLint name));
void pick_handle(GLint n);
static void
zprPick(GLdouble x, GLdouble y,GLdouble delX, GLdouble delY);

void glutResize (int w, int h);


//----------------------------------------------
class coord_
{
public:
  double x,y,z;
  bool flag;

  coord_();
  coord_(double x1, double y1, double z1);
  coord_(const coord_ & c);

  friend ostream& operator<< (ostream &, coord_ &);
  bool empty();

  coord_   operator+(const coord_ &C); 
  coord_   operator*(double a); 

};
//----------------------------------------------

class color_
{
public:
  double x,y,z,a;

  color_ ();
  color_(double x1, double y1, double z1);

  color_(double x1, double y1, double z1, double a1);

  friend ostream& operator<< (ostream &, color_ &);
};
//-------------------------------------------------


class color_map_
{
  map <unsigned int, color_> color_map;
  map <unsigned int, bool> mutable_;
  color_ DEFAULT_COLOR_;

public:
  /*! Constructor */
  color_map_();
  color_map_(map<unsigned int, 
	     map<unsigned int, 
	     double> > cmat);

  color_map_(map<unsigned int, 
	     map<unsigned int, 
	     double> > cmat, 
	     color_ def);

  /*! set alpha */
  void set_alpha(double alpha_factor);

  /*! Accessor */
  color_ &get(unsigned int C);
};
//-------------------------------------------------

class atom_
{
private:

  string a_name;
  coord_ a_coord;
  unsigned int a_chain_pos;
  string aa;
  bool a_visible;
  unsigned int a_attrib;

public:
  /*! Constructors */
  atom_();
  atom_(string atom_name, coord_ atom_coord, unsigned int pos);
  atom_(string atom_name, coord_ atom_coord, unsigned int pos, string aac);
  atom_(string atom_name, coord_ atom_coord, unsigned int pos, unsigned int C);
  atom_(string atom_name, coord_ atom_coord, unsigned int pos, unsigned int C, bool vis);

  atom_(const atom_ &A);

  /*! accessor functions */
  coord_ & coord();
  string name();
  unsigned int pos();
  string amin();
  bool visible();
  unsigned int attrib();

  /*! Function to update atom properties */
  void set_coord(coord_ C);
  void show();
  void hide();
  void set_aa(string aac);
  void set_attrib(unsigned int C);

};
//-------------------------------------------------

class molecule_
{
protected:
  vector <atom_> atoms;
  unsigned int c_length;
  unsigned int n_atoms;

  map <unsigned int, vector <unsigned int> > chain_pos_atoms;
  bool m_visible;

public:
  /*!    Constructor  */
  molecule_();
  molecule_(vector <atom_> A);

  /*!    Accessor  */
  size_t size();
  unsigned int chain_length();
  atom_ & get_atom(unsigned int index);
  bool visible();

/*!    Initialize attribute for atoms  */
  void init_attrib();

/*!    Set attribute for atoms  */
  void set_attrib(unsigned int C);

  void set_attrib(vector <unsigned int> index, unsigned int C);

  void set_attrib_ch(vector <unsigned int> index, unsigned int C);

/*! update atom visibility based on name */
  void hide(string namE);

  void hide_only(vector <string> name_v);

  void show(string namE);

  void show_only(vector <string> name_v);

  /*! update molecule  visibility  */
  void hide();

  /*! update  molecule  visibility  */
  void show();
  
  /*! update atom visibility based on amin */
  void hide_a(string namE);
  /*! update atom visibility based on amin */

  void show_a(string namE);
};
//--------------------------------------------------------
void drawSphere (coord_ p, color_ &color);
void drawStick(coord_ p1,coord_ p2);
void drawStick(coord_ p1,coord_ p2, color_ color);
void drawConnx(coord_ p1, coord_ p2, color_ col);
void drawConnx(coord_ p1, coord_ p2);
void drawAxis (void);
void glutDisplay (void);
void glutIdle (void);
void glutKeyboard (unsigned char key, int x, int y);

/*! Read pdb file */
void data_reader(string pdb_file,
		 string flag_line,
		 vector < molecule_ > & M,
		 string specst="");


//----------------------------------

namespace UTIL_
{
  void proc_spec_str(string,vector<unsigned int>& ATTRIB_CHAIN_INDEX,
			  vector<unsigned int>& ATTRIB_CHAIN_ATTRIB,
			  map<unsigned int,map<unsigned int,unsigned int> >& ATTRIB_CHAIN_POS);
  vector<string> set2vec(set<string>&);

};








