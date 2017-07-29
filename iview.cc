#include "iview.h"


string GL2PS_EXT=".pdf";
GLint GL2PS_TYPE=GL2PS_PDF;

//-------------------------------------------------

const string VERSION="\nSMASH v1.31415 \nCopyright Ishanu Chattopadhyay 2017 UChicago";
const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";
string snapshotfile="";

//-------------------------------------------------
const double r2d = 180/3.14159265359;
bool VERBOSE_=true;
//-------------------------------------------------
//-------------------------------------------------
string APP_NAME="iView v0.1 : ";
/* NB. OpenGL Matrices are COLUMN major. */
vector <string> str_CA;
vector <string> str_CA_CB;
vector <string> str_CA_CB_N;
vector <string> str_CA_CB_N_O;
vector <string> str_ALL;

set <string> ALL_ATOM_NAMES;

double SCALE_=1.0;

double MAX_SEPARATION=5.0;
static string SESSION_ID_="";
static float ALPHA_=1.0;
static GLint IDLE_CALL_TIME_=0;
static int  _mouseX      = 0;
static int  _mouseY      = 0;
static bool _mouseLeft   = false;
static bool _mouseMiddle = false;
static bool _mouseRight  = false;

static double _dragPosX  = 0.0;
static double _dragPosY  = 0.0;
static double _dragPosZ  = 0.0;

static double _matrix[16];
static double _matrixInverse[16];
//-------------------------------------------------
static double _left   = 0.0;
static double _right  = 0.0;
static double _bottom = 0.0;
static double _top    = 0.0;
static double _zNear  = -10.0;
static double _zFar   = 10.0;

static bool SHOW_CONNX=false;
static vector < pair< pair<int, unsigned int>, pair<int, unsigned int> > > CONNX_;

bool BACKGROUND_BLACK=true;
set <double> C_VIOLATION;

//-------------------------------------------------
//-------------------------------------------------
GLfloat ReferencePoint[4] = { 0,0,0,0 };
map <GLint,pair<unsigned int, unsigned int> > NAME_MAP_;
bool update_name_map=true;

//-------------------------------------------------
//-------------------------------------------------
vector<string> set2vec(set<string>& S)
{
  vector<string> vec;
  for(set<string>::iterator itr=S.begin();
      itr!=S.end();
      ++itr)
    vec.push_back(*itr);

  return vec;
}
//-------------------------------------------------
static double
vlen(double x,double y,double z)
{
  return sqrt(x*x+y*y+z*z);
};
//-------------------------------------------------
static void pos(double *px,
		double *py,
		double *pz,
		const int x,
		const int y,
		const int *viewport)
{
  /*
    Use the ortho projection and viewport information
    to map from mouse co-ordinates back into world
    co-ordinates
  */

  *px = (double)(x-viewport[0])/(double)(viewport[2]);
  *py = (double)(y-viewport[1])/(double)(viewport[3]);

  *px = _left + (*px)*(_right-_left);
  *py = _top  + (*py)*(_bottom-_top);
  *pz = _zNear;
};
//-------------------------------------------------
static void getMatrix()
{
  glGetDoublev(GL_MODELVIEW_MATRIX,_matrix);
  invertMatrix(_matrix,_matrixInverse);
};
//-------------------------------------------------
/*
 * From Mesa-2.2\src\glu\project.c
 * Compute the inverse of a 4x4 matrix.  Contributed by scotter@lafn.org
 */
//-------------------------------------------------
static void invertMatrix(const GLdouble *m, 
			 GLdouble *out )
{
#define MAT(m,r,c) (m)[(c)*4+(r)]
#define m11 MAT(m,0,0)
#define m12 MAT(m,0,1)
#define m13 MAT(m,0,2)
#define m14 MAT(m,0,3)
#define m21 MAT(m,1,0)
#define m22 MAT(m,1,1)
#define m23 MAT(m,1,2)
#define m24 MAT(m,1,3)
#define m31 MAT(m,2,0)
#define m32 MAT(m,2,1)
#define m33 MAT(m,2,2)
#define m34 MAT(m,2,3)
#define m41 MAT(m,3,0)
#define m42 MAT(m,3,1)
#define m43 MAT(m,3,2)
#define m44 MAT(m,3,3)
  GLdouble det;
  GLdouble d12, d13, d23, d24, d34, d41;
  GLdouble tmp[16]; /* Allow out == in. */

  /* Inverse = adjoint / det. (See linear algebra texts.)*/
  /* pre-compute 2x2 dets for last two rows when computing */
  /* cofactors of first two rows. */
  d12 = (m31*m42-m41*m32);
  d13 = (m31*m43-m41*m33);
  d23 = (m32*m43-m42*m33);
  d24 = (m32*m44-m42*m34);
  d34 = (m33*m44-m43*m34);
  d41 = (m34*m41-m44*m31);

  tmp[0] =  (m22 * d34 - m23 * d24 + m24 * d23);
  tmp[1] = -(m21 * d34 + m23 * d41 + m24 * d13);
  tmp[2] =  (m21 * d24 + m22 * d41 + m24 * d12);
  tmp[3] = -(m21 * d23 - m22 * d13 + m23 * d12);

  /* Compute determinant as early as possible using these cofactors. */
  det = m11 * tmp[0] + m12 * tmp[1] + m13 * tmp[2] + m14 * tmp[3];
  /* Run singularity test. */
  if (det == 0.0) 
    {
      /* printf("invert_matrix: Warning: Singular matrix.\n"); */
      /*    memcpy(out,_identity,16*sizeof(double)); */
    }
  else 
    {
      GLdouble invDet = 1.0 / det;
      /* Compute rest of inverse. */
      tmp[0] *= invDet;
      tmp[1] *= invDet;
      tmp[2] *= invDet;
      tmp[3] *= invDet;
      
      tmp[4] = -(m12 * d34 - m13 * d24 + m14 * d23) * invDet;
      tmp[5] =  (m11 * d34 + m13 * d41 + m14 * d13) * invDet;
      tmp[6] = -(m11 * d24 + m12 * d41 + m14 * d12) * invDet;
      tmp[7] =  (m11 * d23 - m12 * d13 + m13 * d12) * invDet;
      
      /* Pre-compute 2x2 dets for first two rows when computing */
      /* cofactors of last two rows. */
      d12 = m11*m22-m21*m12;
      d13 = m11*m23-m21*m13;
      d23 = m12*m23-m22*m13;
      d24 = m12*m24-m22*m14;
      d34 = m13*m24-m23*m14;
      d41 = m14*m21-m24*m11;
      
      tmp[8] =  (m42 * d34 - m43 * d24 + m44 * d23) * invDet;
      tmp[9] = -(m41 * d34 + m43 * d41 + m44 * d13) * invDet;
      tmp[10] =  (m41 * d24 + m42 * d41 + m44 * d12) * invDet;
      tmp[11] = -(m41 * d23 - m42 * d13 + m43 * d12) * invDet;
      tmp[12] = -(m32 * d34 - m33 * d24 + m34 * d23) * invDet;
      tmp[13] =  (m31 * d34 + m33 * d41 + m34 * d13) * invDet;
      tmp[14] = -(m31 * d24 + m32 * d41 + m34 * d12) * invDet;
      tmp[15] =  (m31 * d23 - m32 * d13 + m33 * d12) * invDet;
      
      memcpy(out, tmp, 16*sizeof(GLdouble));
    }

#undef m11
#undef m12
#undef m13
#undef m14
#undef m21
#undef m22
#undef m23
#undef m24
#undef m31
#undef m32
#undef m33
#undef m34
#undef m41
#undef m42
#undef m43
#undef m44
#undef MAT
};
//-------------------------------------------------
static void (*selection)(void) = NULL;
static void (*pick)(GLint name) = NULL;
//-------------------------------------------------
void zprSelectionFunc(void (*f)(void))
{
  selection = f;
};
//-------------------------------------------------
void zprPickFunc(void (*f)(GLint name))
{
  pick = f;
};
//-------------------------------------------------
//-------------------------------------------------
/* Draw in selection mode */
static void
zprPick(GLdouble x, GLdouble y,GLdouble delX, GLdouble delY)
{
  GLuint buffer[1024];
  const int bufferSize = sizeof(buffer)/sizeof(GLuint);

  update_name_map=true;
  GLint    viewport[4];
  GLdouble projection[16];

  GLint hits;
  GLint i,j,k;

  GLint  min  = -1;
  GLuint minZ = -1;

  glSelectBuffer(bufferSize,buffer);              /* Selection buffer for hit records */
  glRenderMode(GL_SELECT);                        /* OpenGL selection mode            */
  glInitNames();                                  /* Clear OpenGL name stack          */

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();                                 /* Push current projection matrix   */
  glGetIntegerv(GL_VIEWPORT,viewport);            /* Get the current viewport size    */
  glGetDoublev(GL_PROJECTION_MATRIX,projection);  /* Get the projection matrix        */
  glLoadIdentity();                               /* Reset the projection matrix      */
  gluPickMatrix(x,y,delX,delY,viewport);          /* Set the picking matrix           */
  glMultMatrixd(projection);                      /* Apply projection matrix          */

  glMatrixMode(GL_MODELVIEW);

  if (selection)
    selection();                                 /* Draw the scene in selection mode */

  hits = glRenderMode(GL_RENDER);                 /* Return to normal rendering mode  */

  /* Diagnostic output to stdout */

  if(DEBUG_)
    if (hits!=0)
      {
	printf("hits = %d\n",hits);
	
	for (i=0,j=0; i<hits; i++)
	  {
	    printf("\tsize = %u, min = %u, max = %u : ",buffer[j],buffer[j+1],buffer[j+2]);
	    for (k=0; k < (GLint) buffer[j]; k++)
	      printf("%u ",buffer[j+3+k]);
	    printf("\n");
	    
	    j += 3 + buffer[j];
	  }
      }
  
  /* Determine the nearest hit */
  if (hits)
    {
      for (i=0,j=0; i<hits; i++)
	{
	  if (buffer[j+1]<minZ)
	    {
	      /* If name stack is empty, return -1                */
	      /* If name stack is not empty, return top-most name */

	      if (buffer[j]==0)
		min = -1;
	      else
		min  = buffer[j+2+buffer[j]];

	      minZ = buffer[j+1];
	    }

	  j += buffer[j] + 3;
	}
      pick_handle(min);
    }
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();                         /* Restore projection matrix           */
  glMatrixMode(GL_MODELVIEW);
  update_name_map=false;

  if (pick)
    pick(min);                          /* Pass pick event back to application */
}
//----------------------------------------------
void glutMotion(int x, int y) 
{
  bool changed = false;

  const int dx = x - _mouseX;
  const int dy = y - _mouseY;

  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT,viewport);

  if (dx==0 && dy==0)
    return;

  if ( (_mouseLeft && _mouseRight))
    {
      SCALE_ = exp((double)dy*0.01);

      glTranslatef( ReferencePoint[0], ReferencePoint[1], ReferencePoint[2]);
      glScalef(SCALE_,SCALE_,SCALE_);
      glTranslatef(-ReferencePoint[0],-ReferencePoint[1],-ReferencePoint[2]);
      changed = true;
    }
  else
    if (_mouseLeft)
      {
	double ax,ay,az;
	double bx,by,bz;
	double angle;

	ax = dy;
	ay = dx;
	az = 0.0;
	angle = vlen(ax,ay,az)/(double)(viewport[2]+1)*180.0;

	/* Use inverse matrix to determine local axis of rotation */
	bx = _matrixInverse[0]*ax + _matrixInverse[4]*ay + _matrixInverse[8] *az;
	by = _matrixInverse[1]*ax + _matrixInverse[5]*ay + _matrixInverse[9] *az;
	bz = _matrixInverse[2]*ax + _matrixInverse[6]*ay + _matrixInverse[10]*az;

	glTranslatef( ReferencePoint[0], ReferencePoint[1], ReferencePoint[2]);
	glRotatef(angle,bx,by,bz);
	glTranslatef(-ReferencePoint[0],-ReferencePoint[1],-ReferencePoint[2]);
	changed = true;
      }
    else
      if (_mouseRight)
	{
	  double px,py,pz;
	  pos(&px,&py,&pz,x,y,viewport);
	  glLoadIdentity();
	  glTranslatef(px-_dragPosX,py-_dragPosY,pz-_dragPosZ);
	  glMultMatrixd(_matrix);

	  _dragPosX = px;
	  _dragPosY = py;
	  _dragPosZ = pz;
	  changed = true;
	}
  _mouseX = x;
  _mouseY = y;

  if (changed)
    {
      getMatrix();
      glutPostRedisplay();
    }
}
//----------------------------------------------------------------
void glutMouse(int button, int state, int x, int y) 
{
  GLint viewport[4];

  /* Do picking */
  if ((state==GLUT_DOWN) && (button == GLUT_MIDDLE_BUTTON))
    zprPick(x,glutGet(GLUT_WINDOW_HEIGHT)-1-y,3,3);

  _mouseX = x;
  _mouseY = y;

  if (state==GLUT_UP)
    switch (button)
      {
      case GLUT_LEFT_BUTTON:   _mouseLeft   = false; break;
      case GLUT_MIDDLE_BUTTON: _mouseMiddle = false; break;
      case GLUT_RIGHT_BUTTON:  _mouseRight  = false; break;
      }
  else
    switch (button)
      {
      case GLUT_LEFT_BUTTON:   _mouseLeft   = true; break;
      case GLUT_MIDDLE_BUTTON: _mouseMiddle = true; break;
      case GLUT_RIGHT_BUTTON:  _mouseRight  = true; break;
      }

  glGetIntegerv(GL_VIEWPORT,viewport);
  pos(&_dragPosX,&_dragPosY,&_dragPosZ,x,y,viewport);
  glutPostRedisplay();
}
//----------------------------------------------
void glutResize (int w, int h)
{	
  if (!h)
    return;
  glViewport(0,0,w,h);

  _top    =  1.0;
  _bottom = -1.0;
  _left   = -(double)w/(double)h;
  _right  = -_left;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(_left,_right,_bottom,_top,_zNear,_zFar);

  glMatrixMode(GL_MODELVIEW);
  glutPostRedisplay();
}
//----------------------------------------------
//----------------------------------------------
coord_::coord_()
{
  flag=false;
};

coord_::coord_(double x1, double y1, double z1)
{
  x=x1; y=y1; z=z1;flag=true;
}

coord_::coord_(const coord_ & c)
{
  x=c.x;y=c.y;z=c.z;flag=true;
};

bool coord_::empty()
{
  return flag;
}

coord_   coord_::operator+(const coord_ &C)
{
  coord_  c(x+C.x, y+C.y,z+C.z);
  return c;
}; 

coord_   coord_::operator*(double a)
{
  coord_  c(x*a,y*a,z*a);
  return c;
}; 

//-------------------------------------------------
ostream& operator << (ostream &out, coord_ &C)
{
  out << C.x << " " << C.y << " " << C.z;
  return out;
};
//-------------------------------------------------
ostream& operator << (ostream &out, vector<string> &C)
{
  for(vector<string>::iterator itr=C.begin();
      itr!=C.end();
      ++itr)
    out << *itr << " ";
  return out;
};
//-------------------------------------------------
color_::color_()
{};

color_::color_(double x1, double y1, double z1)
{ 
  x=x1;
  if (x1 > 1.0)
    x=1.0;
  if (x1 < 0)
    x=0.0;

  y=y1;
  if (y1 > 1.0)
    y=1.0;
  if (y1 < 0)
    y=0.0;

  z=z1;
  if (z1 > 1.0)
    z=1.0;
  if (z1 < 0)
    z=0.0;

  a=ALPHA_;
  if (a > 1.0)
    a = 1.0;
};

color_::color_(double x1, double y1, double z1, double a1)
{ 
  x=x1;
  if (x1 > 1.0)
    x=1.0;
  if (x1 < 0)
    x=0.0;

  y=y1;
  if (y1 > 1.0)
    y=1.0;
  if (y1 < 0)
    y=0.0;

  z=z1;
  if (z1 > 1.0)
    z=1.0;
  if (z1 < 0)
    z=0.0;

  a=a1;
  if (a > 1.0)
    a = 1.0;
};

//-------------------------------------------------
ostream& operator<< (ostream & out, color_ & col)
{
  out << col.x << " " << col.y << " " << col.z;
  return out;
};
//-------------------------------------------------

color_map_::color_map_(){};
color_map_::color_map_(map<unsigned int, 
		       map<unsigned int, 
		       double> > cmat)
{
  for (map<unsigned int, 
	 map<unsigned int,
	 double> >::iterator itr=cmat.begin();
       itr != cmat.end();
       ++itr)
    if ((itr->second.find(0) != itr->second.end())
	&& (itr->second.find(1) != itr->second.end())
	&& (itr->second.find(2) != itr->second.end())
	)
      {
	if ((itr->second.find(3) == itr->second.end()))
	  color_map[itr->first] = 
	    color_ ((itr->second)[0],
		    (itr->second)[1],
		    (itr->second)[2]); 
	else
	  color_map[itr->first] = 
	    color_ ((itr->second)[0],
		    (itr->second)[1],
		    (itr->second)[2], 
		    (itr->second)[3]); 
      }

  for (map<unsigned int, 
	 map<unsigned int, 
	 double> >::iterator itr=cmat.begin();
       itr != cmat.end();
       ++itr)
    if ((itr->second.find(4) != itr->second.end()))
      mutable_[itr->first] = (itr->second)[4];
    else
      mutable_[itr->first] = false;

  DEFAULT_COLOR_ = color_ (.2,.2,.2);
};

color_map_::color_map_(map<unsigned int, 
		       map<unsigned int, 
		       double> > cmat, 
		       color_ def)
{
  for (map<unsigned int,
	 map<unsigned int,
	 double> >::iterator itr=cmat.begin();
       itr != cmat.end();
       ++itr)
    if ((itr->second.find(0) != itr->second.end())
	&& (itr->second.find(1) != itr->second.end())
	&& (itr->second.find(2) != itr->second.end())
	)
      {
	if ((itr->second.find(3) == itr->second.end()))
	  color_map[itr->first] = 
	    color_ ((itr->second)[0],
		    (itr->second)[1],
		    (itr->second)[2]); 
	else
	  color_map[itr->first] = 
	    color_ ((itr->second)[0],
		    (itr->second)[1],
		    (itr->second)[2], 
		    (itr->second)[3]); 
      }

  for (map<unsigned int, 
	 map<unsigned int, 
	 double> >::iterator itr=cmat.begin();
       itr != cmat.end();
       ++itr)
    if ((itr->second.find(4) != itr->second.end()))
      mutable_[itr->first] = (itr->second)[4];
    else
      mutable_[itr->first] = false;

  DEFAULT_COLOR_ = def;
};

/*! set alpha */
void color_map_::set_alpha(double alpha_factor)
{
  for (map <unsigned int, color_>::iterator itr=color_map.begin();
       itr!=color_map.end();
       ++itr)
    if (mutable_[itr->first])
      {
	if (itr->second.a * alpha_factor < 1.0)
	  itr->second.a *= alpha_factor;
	else 
	  itr->second.a = 1.0;
      }
};

/*! Accessor */
color_ & color_map_::get(unsigned int C)
{
  if (color_map.find(C) != color_map.end())
    return color_map[C];
  else
    return DEFAULT_COLOR_;
};


//-------------------------------------------------
color_map_ *COLOR;
//-------------------------------------------------

string IMGFILE;
unsigned int IMG_COUNT=0, HR_IMG_COUNT=0;

GLubyte * image_buffer; //RGB bits  

/* window width and height */
int winW = 640;
int winH = 480;

/* old position of the mouse */
int oldX = -13;
int oldY = -13;

/* mouse state, UP or DOWN */
int mState = UP;

/* current axis of rotation */
int axisRot = X_AXIS;

/* amount to rotate about axis */
float rotate_ = 0.0;

/* vector which describes the axis to rotate about */
vec3_t axis = {0.0, 1.0, 0.0};

/* global rotation, for use with the mouse */
vec3_t gRot = {0,0,0};
float rotspeed = 0.005f;

/* should we animate ? */
int aniOn = 0;

double zoom_ = 1.0;
double r1 = 0.2;
double r2 = 0.2;
double SPHERE_RAD_ = 1;
int SPHERE_RES_ = 30;
int precision = SPHERE_RES_;

string configfile;

//-------------------------------------------------
//-------------------------------------------------

const GLfloat light_ambient[]  = { 0.0f, 0.0f, 0.0f, 1.0f };
const GLfloat light_diffuse[]  = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_position[] = { 2.0f, 5.0f, 5.0f, 0.0f }; 
 
const GLfloat mat_ambient[]    = { 0.7f, 0.7f, 0.7f, 1.0f };
const GLfloat mat_diffuse[]    = { 0.8f, 0.8f, 0.8f, 1.0f };
const GLfloat mat_specular[]   = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat high_shininess[] = { 100.0f }; 

//-------------------------------------------------
//-------------------------------------------------
// atoms are defined by chain pos, attrib flag, 

atom_::atom_(){};
atom_::atom_(string atom_name, coord_ atom_coord, unsigned int pos)
{
  a_name = atom_name;
  a_coord = coord_ (atom_coord.x, atom_coord.y,atom_coord.z);
  a_chain_pos = pos;
  a_visible = true;
  aa = "";
  a_attrib = 0;
};
atom_::atom_(string atom_name, coord_ atom_coord, unsigned int pos, string aac)
{
  a_name = atom_name;
  a_coord = coord_ (atom_coord.x, atom_coord.y,atom_coord.z);
  a_chain_pos = pos;
  aa = aac;
  a_visible = true;
  a_attrib = 0;
};
atom_::atom_(string atom_name, coord_ atom_coord, unsigned int pos, unsigned int C)
{
  a_name = atom_name;
  a_coord = coord_ (atom_coord.x, atom_coord.y,atom_coord.z);
  a_chain_pos = pos;
  a_visible = true;
  aa = "";
  a_attrib = C;
};
atom_::atom_(string atom_name, coord_ atom_coord, unsigned int pos, unsigned int C, bool vis)
{
  a_name = atom_name;
  a_coord = coord_ (atom_coord.x, atom_coord.y,atom_coord.z);
  a_chain_pos = pos;
  a_visible = vis;
  aa = "";
  a_attrib = C;
};

atom_::atom_(const atom_ &A)
{
  a_name = A.a_name;
  a_coord = A.a_coord;
  a_chain_pos = A.a_chain_pos;
  aa = A.aa;
  a_visible = A.a_visible;
  a_attrib = A.a_attrib;
};

/*! accessor functions */
coord_ & atom_::coord(){return a_coord;};
string atom_::name(){return a_name;};
unsigned int atom_::pos(){return a_chain_pos;};
string atom_::amin(){return aa;};
bool atom_::visible(){return a_visible;};
unsigned int atom_::attrib(){return a_attrib;};

/*! Function to update atom properties */
void atom_::set_coord(coord_ C){a_coord = C;};
void atom_::show(){a_visible = true;};
void atom_::hide(){a_visible = false;};
void atom_::set_aa(string aac){aa = aac;};
void atom_::set_attrib(unsigned int C){a_attrib = C;};

//--------------------------------------------------------------
//------------------------------------------

molecule_::molecule_(){};
molecule_::molecule_(vector <atom_> A)
{
  m_visible=true;
  atoms = A;
  n_atoms = atoms.size();
  set <unsigned int> cpos;
  for (vector <atom_>::iterator itr=atoms.begin();
       itr!= atoms.end();
       ++itr)
    cpos.insert(itr->pos());
  c_length=cpos.size();

  for (unsigned int i=0; i < n_atoms;++i)
    chain_pos_atoms[atoms[i].pos()].push_back(i);
};

/*!    Accessor  */
size_t molecule_::size(){return n_atoms;};
unsigned int molecule_::chain_length(){return c_length;};
atom_ & molecule_::get_atom(unsigned int index)
{
  if (index > n_atoms)
    index = n_atoms;
  return atoms[index];
}
bool molecule_::visible(){return m_visible;};

/*!    Initialize attribute for atoms  */
void molecule_::init_attrib()
{
  for(vector <atom_>::iterator itr=atoms.begin();
      itr != atoms.end();
      ++itr)
    itr->set_attrib(0);
};

/*!    Initialize attribute for atoms  */
void molecule_::set_attrib(unsigned int C)
{
  for(vector <atom_>::iterator itr=atoms.begin();
      itr != atoms.end();
      ++itr)
    itr->set_attrib(C);
};

/*!    Set attribute for atoms  */
void molecule_::set_attrib(vector <unsigned int> index, unsigned int C)
{
  for (unsigned int i=0; i < index.size();++i)
    if (index[i] < n_atoms)
      atoms[index[i]].set_attrib(C);
};

/*! Set attribute for atoms  */
void molecule_::set_attrib_ch(vector <unsigned int> index, unsigned int C)
{

  for (unsigned int i=0; i < index.size();++i)
    for (unsigned int j=0; j < chain_pos_atoms[index[i]].size(); ++j)
      atoms[chain_pos_atoms[index[i]][j]].set_attrib(C);
};

/*! update atom visibility based on name */
void molecule_::hide(string namE)
{
  for(vector <atom_>::iterator itr=atoms.begin();
      itr != atoms.end();
      ++itr)
    if (itr->name() == namE)
      itr->hide();
};

/*! update atom visibility based on name */
void molecule_::hide_only(vector <string> name_v)
{
  set <string> s_name_v(name_v.begin(), name_v.end());
  for(vector <atom_>::iterator itr=atoms.begin();
      itr != atoms.end();
      ++itr)
    if (s_name_v.find(itr->name()) != s_name_v.end() )
      itr->hide();
    else
      itr->show();
};

/*! update atom visibility based on name */
void molecule_::show(string namE)
{
  for(vector <atom_>::iterator itr=atoms.begin();
      itr != atoms.end();
      ++itr)
    if (itr->name() == namE)
      itr->show();
};

/*! update atom visibility based on name */
void molecule_::show_only(vector <string> name_v)
{
  set <string> s_name_v(name_v.begin(), name_v.end());
  for(vector <atom_>::iterator itr=atoms.begin();
      itr != atoms.end();
      ++itr)
    if (s_name_v.find(itr->name()) != s_name_v.end() )
      itr->show();
    else
      itr->hide();
};

/*! update molecule  visibility  */
void molecule_::hide(){m_visible=false;  };

/*! update  molecule  visibility  */
void molecule_::show(){m_visible=true;  };
  
/*! update atom visibility based on amin */
void molecule_::hide_a(string namE)
{
  for(vector <atom_>::iterator itr=atoms.begin();
      itr != atoms.end();
      ++itr)
    if (itr->amin() == namE)
      itr->hide();
};
/*! update atom visibility based on amin */
void molecule_::show_a(string namE)
{
  for(vector <atom_>::iterator itr=atoms.begin();
      itr != atoms.end();
      ++itr)
    if (itr->amin() == namE)
      itr->show();
};



//--------------------------------------------------------
vector <molecule_ > M;
//-------------------------------------------------------
void drawSphere (coord_ p, color_ &color)
{
  glPushMatrix ();
  glScalef(zoom_, zoom_, zoom_);
  glTranslatef (p.x,p.y,p.z);
  glColor4f (color.x,color.y,color.z, color.a);
  glutSolidSphere (SPHERE_RAD_, SPHERE_RES_, SPHERE_RES_);
  glPopMatrix ();
}
//-------------------------------------------------
void drawStick(coord_ p1,coord_ p2)
{
  double vx=p2.x-p1.x;
  double vy=p2.y-p1.y;
  double vz=p2.z-p1.z;
  double c_height=sqrt(vx*vx+vy*vy+vz*vz);


  if (c_height>MAX_SEPARATION) 
    {
      if (C_VIOLATION.insert(c_height).second)
	cout << "MAX C VIOLATION WARNING " << c_height << endl;
      return;
    }

  double  rx=-vy*vz,  ry=vx*vz, ax=0.0;
  if (vz==0)
    {
      ax =-r2d*acos(vx/c_height);  // rotation angle in x-y plane
      if (vx<=0)
	ax=-ax;
    }    
  else
    {
      ax=r2d*acos(vz/c_height);	// rotation angle
      if (vz<=0)
	ax=-ax;
    }
  glPushMatrix ();
  glScalef(zoom_, zoom_, zoom_);
  glTranslatef(p1.x,p1.y,p1.z);
  if (vz==0)
    {
      glRotated(90.0, 0, 1, 0.0);		// Rotate & align with x axis
      glRotated(ax, -1.0, 0.0, 0.0);		// Rotate to point 2 in x-y plane
    }
  else 
    glRotated(ax, rx, ry, 0);			// Rotate about rotation vector
  GLUquadricObj *quadObj = gluNewQuadric();
  gluQuadricNormals(quadObj, GLU_SMOOTH);
  gluCylinder(quadObj,r1,r2,c_height,precision,precision);	// Draw A cylinder
  glPopMatrix ();
  gluDeleteQuadric(quadObj);
}

//-------------------------------------------------
void drawStick(coord_ p1,coord_ p2, color_ color)
{
  double vx=p2.x-p1.x;
  double vy=p2.y-p1.y;
  double vz=p2.z-p1.z;
  double c_height=sqrt(vx*vx+vy*vy+vz*vz);

  if (c_height>MAX_SEPARATION) 
    {
      if (C_VIOLATION.insert(c_height).second)
	cout << "MAX C VIOLATION WARNING " << c_height << endl;
      return;
    }

  double  rx=-vy*vz,  ry=vx*vz, ax=0.0;
  if (vz==0)
    {
      ax =-r2d*acos(vx/c_height);  // rotation angle in x-y plane
      if (vx<=0)
	ax=-ax;
    }    
  else
    {
      ax=r2d*acos(vz/c_height);	// rotation angle
      if (vz<=0)
	ax=-ax;
    }
  glPushMatrix ();
  glColor3f (color.x,color.y,color.z);
  glScalef(zoom_, zoom_, zoom_);
  glTranslatef(p1.x,p1.y,p1.z);
  if (vz==0)
    {
      glRotated(90.0, 0, 1, 0.0);	 // Rotate & align with x axis
      glRotated(ax, -1.0, 0.0, 0.0);	 // Rotate to point 2 in x-y plane
    }
  else 
    glRotated(ax, rx, ry, 0);		 // Rotate about rotation vector
  GLUquadricObj *quadObj = gluNewQuadric();
  gluQuadricNormals(quadObj, GLU_SMOOTH);
  gluCylinder(quadObj,r1,r2,c_height,precision,precision);   // Draw A cylinder
  glPopMatrix ();
  gluDeleteQuadric(quadObj);
}
//-------------------------------------------------
void drawConnx(coord_ p1, coord_ p2, color_ col)
{
  glBegin (GL_LINES);
  glColor3f (col.x,col.y,col.z);
  glVertex3f (p1.x,p1.y,p1.z);
  glVertex3f (p2.x,p2.y,p2.z);
  glEnd ();
}
//-------------------------------------------------
void drawConnx(coord_ p1, coord_ p2)
{
  glBegin (GL_LINES);
  glColor3f (COLOR->get(1000000).x,COLOR->get(1000000).y,COLOR->get(1000000).z);
  glVertex3f (p1.x,p1.y,p1.z);
  glVertex3f (p2.x,p2.y,p2.z);
  glEnd ();
}
//-------------------------------------------------
/*   drawAxis */
//-------------------------------------------------
void drawAxis (void)
{
  glColor3f (0.5, 0.5, 0.5);
  glBegin (GL_LINES);
  glColor3f (0.5, 0.0, 0.0);
  glVertex3f (-20.0, 0.0, 0.0);
  glVertex3f (20.0, 0.0, 0.0);

  glColor3f (0.0, 0.5, 0.0);
  glVertex3f (0.0, 20.0, 0.0);
  glVertex3f (0.0, -20.0, 0.0);

  glColor3f (0.0, 0.0, 0.5);
  glVertex3f (0.0, 0.0, -20.0);
  glVertex3f (0.0, 0.0, 20.0);
  glEnd ();
}
//-------------------------------------------------
/*  glutDisplay*/ 
void glutDisplay (void)
{
  GLint curr_time=glutGet(GLUT_ELAPSED_TIME);
  if (aniOn)
    {     
      glTranslatef( ReferencePoint[0], ReferencePoint[1], ReferencePoint[2]);
      glRotatef (rotspeed*(curr_time-IDLE_CALL_TIME_), axis[0], axis[1], axis[2]);
      glTranslatef(-ReferencePoint[0],-ReferencePoint[1],-ReferencePoint[2]);
      getMatrix();
      glutPostRedisplay ();
    }
  IDLE_CALL_TIME_=curr_time;
  GLint objcount=0;
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  for (unsigned int i=0; i < M.size(); ++i)
    if (M[i].visible())
      {
	bool flag=false;
	coord_ last_atom_coord;
	for (unsigned int j=0; j < M[i].size(); ++j)
	  if ( M[i].get_atom(j).visible()  )
	    {
	      if (update_name_map)
		{
		  NAME_MAP_[objcount] = make_pair(i,j);
		  glPushName(objcount++);
		}
	      drawSphere (M[i].get_atom(j).coord(),
			  COLOR->get(M[i].get_atom(j).attrib()) );
	      if (update_name_map)
		glPopName();
	      if (flag)
		drawStick(M[i].get_atom(j).coord(),last_atom_coord);
	      flag=true;
	      last_atom_coord = M[i].get_atom(j).coord();
	    }
      }
  if (SHOW_CONNX)
    for (vector < pair< pair<int, unsigned int>, 
	   pair<int, unsigned int> > >::iterator itr=CONNX_.begin();
	 itr != CONNX_.end();
	 ++itr)
      drawConnx(M[itr->first.first].get_atom(itr->first.second).coord(),
		M[itr->second.first].get_atom(itr->second.second).coord());
  /*glPushMatrix ();
    gleSetJoinStyle (TUBE_NORM_EDGE | TUBE_JN_ANGLE | TUBE_JN_CAP);
    gleHelicoid (1.0, 6.0, 0.0, -3.0, 4.0, 0x0, 0x0, 0.0, 1080.0);
    glPopMatrix ();
  */
  glutSwapBuffers();
  // glutPostRedisplay();
  //  update_name_map=false;
}
//-------------------------------------------------------
void pick_handle(GLint n)
{
  map <GLint,pair<unsigned int, unsigned int> >::iterator itr=NAME_MAP_.find(n);

  cout << n << " "
       << itr->second.first
       << " "
       <<  M[itr->second.first].get_atom(itr->second.second).pos()
       << endl;
  // M[NAME_MAP_[n].first].get_atom(NAME_MAP_[n].second).hide();
  //glutPostRedisplay ();
}
//-------------------------------------------------
/*    glutIdle   Idle function, makes it animate. */
void glutIdle (void)
{	
  GLint curr_time=glutGet(GLUT_ELAPSED_TIME);
  if (aniOn)
    {     
      glTranslatef( ReferencePoint[0], ReferencePoint[1], ReferencePoint[2]);
      glRotatef (rotspeed*(curr_time-IDLE_CALL_TIME_), axis[0], axis[1], axis[2]);
      glTranslatef(-ReferencePoint[0],-ReferencePoint[1],-ReferencePoint[2]);
      getMatrix();
      glutPostRedisplay ();
    }
  IDLE_CALL_TIME_=curr_time;
}
//-------------------------------------------------
/*    glutKeyboard   Keyboard handler.*/
void glutKeyboard (unsigned char key, int x, int y)
{
  (void) x;
  (void) y;
  switch (key)
    {
    case 'q':
    case 'Q':
      exit (1);
      break;
    case 'x':
    case 'X':
      /* axis of rotation is now X */
      axisRot = X_AXIS;
      axis[0] = 1.0;
      axis[1] = axis[2] = 0.0;
      break;
    case 'y':
    case 'Y':
      /* axis of rotation is now Y */
      axisRot = Y_AXIS;
      axis[1] = 1.0;
      axis[0] = axis[2] = 0.0;
      break;
    case 'z':
    case 'Z':
      /* axis of rotation is now Z */
      axisRot = Z_AXIS;
      axis[2] = 1.0;
      axis[1] = axis[0] = 0.0;
      break;
    case 'u':
    case 'U':
      /* increase rotspeed */
      rotspeed *= 1.01;
      break;
    case 'd':
    case 'D':
      /* decrease rotspeed */
      rotspeed *= 0.99;
      break;
    case 's':
      /* capture viewport */
      if (glutGetModifiers() != GLUT_ACTIVE_ALT)
	CaptureViewPort();
      else
	{
	  CaptureViewPort_HR();
	}
      break;
    case 't':
      if (glutGetModifiers() != GLUT_ACTIVE_ALT)
	COLOR->set_alpha(1.05);
      else
	COLOR->set_alpha(0.95);
      break;
    case 'c':
      if (glutGetModifiers() == GLUT_ACTIVE_ALT)
	update_show(str_CA);
      break;
    case 'b':
      if (glutGetModifiers() == GLUT_ACTIVE_ALT)
	update_show(str_CA_CB);
      break;
    case 'n':
      if (glutGetModifiers() == GLUT_ACTIVE_ALT)
	update_show(str_CA_CB_N);
      break;
    case 'o':
      if (glutGetModifiers() == GLUT_ACTIVE_ALT)
	update_show(str_CA_CB_N_O);
      break;
    case 'a':
      if (glutGetModifiers() == GLUT_ACTIVE_ALT)
	update_show(str_ALL);
      break;
    case ' ':
      /* animation on/off */
      aniOn = !aniOn;
      break;
    }
  glutPostRedisplay ();
}
//-------------------------------------------------
/*! Read pdb file */
void data_reader(string pdb_file,
		 string flag_line,
		 vector < molecule_ > & M,
		 string SPECSTR)
{
  ifstream DATA(pdb_file.c_str());
  string line;
  vector <atom_ > A;
  unsigned int last_pos=0;

  while (getline (DATA, line))
    {
      string atom_name;
      coord_ atom_coord;
      unsigned int chain_pos;
      string aa;

      if (line.substr(0,4) == flag_line)
	{
	  aa = line.substr(17,3);
	  atom_name =  line.substr(13,2);
	  trim(atom_name);
	  atom_coord = coord_ (atof(line.substr(30,8).c_str()),
			       atof(line.substr(38,8).c_str()),
			       atof(line.substr(46,8).c_str()));
	  chain_pos = atoi(line.substr(22,6).c_str());

	  if ((chain_pos < last_pos ))
	    if (!A.empty())
	      {
		if (DEBUG_)
		  cout << "## ATOMS READ: " << A.size() << " chain_len: " <<  last_pos<< endl;
		M.push_back(molecule_ (A));
		A.clear();
	      }
	  ALL_ATOM_NAMES.insert(atom_name);
	  A.push_back(atom_ (atom_name,atom_coord,chain_pos,aa) );  
	  last_pos = chain_pos;
	}
    }
  if (DEBUG_)
    cout << "## ATOMS READ: " << A.size() << " chain_len: " <<  last_pos << endl;

  M.push_back(molecule_ (A) );

  coord_ S(0,0,0);
  unsigned int Ms=0;
    
  for (unsigned int m=0; m < M.size(); ++m)
    {
      for (unsigned int i=0; i < M[m].size();++i)
	S = S + M[m].get_atom(i).coord();
      Ms += M[m].size();
    }
     
  S = S*(1/(Ms + 0.0));
  coord_ Sm(-S.x,-S.y,-S.z);
  
  for (unsigned int m=0; m < M.size(); ++m)
    for (unsigned int i=0; i < M[m].size();++i)
      M[m].get_atom(i).set_coord(M[m].get_atom(i).coord() + Sm);

  molviewConfig(SPECSTR);


  for (unsigned int i=0; i<M.size();++i)
    cout << "MOLECULE " << i << endl;
};
//--------------------------------------------------------------
/*! Write JPEG file from sceenshot  using libjpeg */
void write_JPEG_file (const char *filename, 
		      int quality,
		      int image_width, 
		      int image_height)
{
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  FILE * outfile;		/* target file */
  JSAMPROW row_pointer[1];	/* pointer to JSAMPLE row[s] */
  int row_stride;		/* physical row width in image buffer */
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);

  if ((outfile = fopen(filename, "wb")) == NULL) 
    {
      fprintf(stderr, "can't open %s\n", filename);
      exit(1);
    }
  jpeg_stdio_dest(&cinfo, outfile);
  cinfo.image_width = image_width; /* image width and height, in pixels */
  cinfo.image_height = image_height;
  cinfo.input_components = 3; /* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB;/* colorspace of input image */
  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo, quality, TRUE); /* limit to baseline-JPEG */
  jpeg_start_compress(&cinfo, TRUE);
  row_stride = image_width * 3;	/* JSAMPLEs per row in image_buffer */

  while (cinfo.next_scanline < cinfo.image_height) 
    {
      row_pointer[0] = 
	& image_buffer[(cinfo.image_height
			-1-cinfo.next_scanline)*row_stride];
      (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }
  jpeg_finish_compress(&cinfo);
  fclose(outfile);
  jpeg_destroy_compress(&cinfo);
};
//-------------------------------------------------
/*! Function to print screen */
void CaptureViewPort()  
{  
  ostringstream img_count;   // stream used for the conversion
  img_count << IMG_COUNT++;    

  int NAMEWIDTH=10;
  string zpad_="";
  for (unsigned int i=0; i<=NAMEWIDTH-img_count.str().length();++i)
    zpad_='0'+zpad_;

  GLint viewport[4];  //get current viewport  
  glGetIntegerv(GL_VIEWPORT, viewport);  
  glFinish(); //finish all commands of OpenGL  
  winW = viewport[2];  
  winH = viewport[3];  
  image_buffer = new GLubyte[winW*3*winH];
  //read pixel from frame buffer  
  glPixelStorei(GL_PACK_ALIGNMENT,1);  
  glPixelStorei(GL_PACK_ROW_LENGTH, 0);  
  glPixelStorei(GL_PACK_SKIP_ROWS, 0);  
  glPixelStorei(GL_PACK_SKIP_PIXELS, 0); 
  
  glReadPixels(0, 0, winW, winH, 
	       GL_RGB, 
	       GL_UNSIGNED_BYTE, 
	       image_buffer);  
  write_JPEG_file ((IMGFILE + 
		    SESSION_ID_+ zpad_+
		    img_count.str() + 
		    ".jpeg").c_str(), 
		   100, winW, winH);
  delete [] image_buffer;
};
//-------------------------------------------------
void CaptureViewPort_HR()  
{
  ostringstream img_count;   // stream used for the conversion
  img_count << HR_IMG_COUNT++;    
  FILE *fp;
  int state = GL2PS_OVERFLOW, buffsize = 1024*4096*10;
  fp = fopen((IMGFILE + 
	      SESSION_ID_+ 
	      img_count.str() + 
	      GL2PS_EXT).c_str(), "wb");
  // printf("Writing 'out.eps'... ");
  while(state == GL2PS_OVERFLOW)
    {
      buffsize += 1024*1024*5;
      gl2psBeginPage("image", 
		     "iView", 
		     NULL, 
		     GL2PS_TYPE, 
		     GL2PS_SIMPLE_SORT,
		     GL2PS_DRAW_BACKGROUND 
		     | GL2PS_USE_CURRENT_VIEWPORT | GL2PS_OCCLUSION_CULL,
		     // | GL2PS_OCCLUSION_CULL | GL2PS_COMPRESS,
		     GL_RGBA, 0, NULL, 0, 0, 0, 
		     buffsize, fp, "out.pdf");
      glutDisplay();
      state = gl2psEndPage();
    }
  fclose(fp);
  printf("Done!\n");
}
//-------------------------------------------------
/*!   glInit   Sets up various OpenGL stuff. */
void glInit() {
  IMGFILE = snapshotfile;
  glClearColor(1,1,1,1);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK); 
 
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS); 
  
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glEnable(GL_LIGHT0);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_LIGHTING); 
 
  glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position); 
 
  glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
  glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
  glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess); 

  if (BACKGROUND_BLACK)
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Set background color to black and opaque
  else
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f); 
  glClearDepth(1.0f);                   // Set background depth to farthest
  // glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
  glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
  glShadeModel(GL_SMOOTH);   // Enable smooth shading
  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Nice perspective corrections
  glScalef(.01,.01,.01);

  str_CA.push_back("CA");;
  str_CA_CB=str_CA;
  str_CA_CB.push_back("C");
  str_CA_CB_N=str_CA_CB;
  str_CA_CB_N.push_back("N");
  str_CA_CB_N_O=str_CA_CB_N;
  str_CA_CB_N_O.push_back("O");
  str_ALL=set2vec(ALL_ATOM_NAMES);

  
  struct timeval tim;  
  gettimeofday(&tim, NULL);  
  // double t0=tim.tv_sec;
  double t1=(tim.tv_usec);  
  ostringstream session_id_;   // stream used for the conversion
  // session_id_ <<  t0;
  //session_id_ << ":";
  session_id_ << t1;
  SESSION_ID_ = "_" + session_id_.str()+ "_";
  if (VERBOSE_)
    cout << SESSION_ID_ << endl;


}
//-------------------------------------------------

void molviewConfig(string specstr)
{
  unsigned int GLOBAL_DEFAULT_ATTRIB;
  vector <unsigned int> DEFAULT_ATTRIB_INDEX;
  vector <unsigned int> DEFAULT_ATTRIB;
  vector <unsigned int> SHOW_ATOMS_INDEX;
  map <unsigned int, map <unsigned int, string > > ATOM_NAMES;
  unsigned int GLOBAL_ATTRIB_CHAIN_ATTRIB = 0;
  vector <unsigned int> GLOBAL_ATTRIB_CHAIN_POS;
  vector <unsigned int> ATTRIB_CHAIN_INDEX;
  vector <unsigned int> ATTRIB_CHAIN_ATTRIB;
  map <unsigned int, map <unsigned int, unsigned int > > ATTRIB_CHAIN_POS;
  vector <unsigned int> VISIBLE_INDEX;
  vector <bool> VISIBLE_FLAG;
  vector <double>  DEFAULT_COLOR(4,0.2);
  map <unsigned int, map <unsigned int, double> > COLOR_MAP;
  bool GLOBAL_VISIBLE_=true;
  string gl2ps_type="GL2PS_PDF";
  string APP_NAME_EXT="";


  CONFIG cfg(configfile);

  cfg.set(SPHERE_RAD_,"SPHERE_RAD");
  cfg.set(SPHERE_RES_,"SPHERE_RES");

  cfg.set(GLOBAL_DEFAULT_ATTRIB,"GLOBAL_DEFAULT_ATTRIB");
  cfg.set_vector<unsigned int>(DEFAULT_ATTRIB_INDEX,"DEFAULT_ATTRIB_INDEX");
  cfg.set_vector<unsigned int>(DEFAULT_ATTRIB,"DEFAULT_ATTRIB");
  cfg.set_vector<unsigned int>(SHOW_ATOMS_INDEX,"SHOW_ATOMS_INDEX");
  cfg.set_map<unsigned int,unsigned int,string>(ATOM_NAMES,"SHOW_ATOMS");
  cfg.set(GLOBAL_ATTRIB_CHAIN_ATTRIB,"GLOBAL_ATTRIB_CHAIN_ATTRIB");
  cfg.set_vector<unsigned int>(GLOBAL_ATTRIB_CHAIN_POS,"GLOBAL_ATTRIB_CHAIN_POS");
  cfg.set_vector<unsigned int>(ATTRIB_CHAIN_INDEX,"ATTRIB_CHAIN_INDEX");
  cfg.set_vector<unsigned int>(ATTRIB_CHAIN_ATTRIB,"ATTRIB_CHAIN_ATTRIB");
  cfg.set_map<unsigned int,unsigned int,
	      unsigned int>(ATTRIB_CHAIN_POS,"ATTRIB_CHAIN_POS");
  cfg.set(GLOBAL_VISIBLE_,"GLOBAL_VISIBLE");
  cfg.set_vector<unsigned int>(VISIBLE_INDEX,"VISIBLE_INDEX");
  cfg.set_vector<bool>(VISIBLE_FLAG,"VISIBLE");
  cfg.set_vector<double>(DEFAULT_COLOR,"DEFAULT_COLOR");
  cfg.set_map<unsigned int,unsigned int,double>(COLOR_MAP,"COLOR_MAP");
  cfg.set(gl2ps_type,"GL2PS_TYPE");
  cfg.set(APP_NAME_EXT,"APP_NAME_EXT");
  cfg.set(BACKGROUND_BLACK,"BACKGROUND_BLACK");

  APP_NAME += APP_NAME_EXT;


  UTIL_::proc_spec_str(specstr,
		       ATTRIB_CHAIN_INDEX,
		       ATTRIB_CHAIN_ATTRIB,
		       ATTRIB_CHAIN_POS);
    


  if (gl2ps_type=="")
    {
      GL2PS_TYPE=GL2PS_PDF;
      GL2PS_EXT=".pdf";
    }
  if (gl2ps_type=="PDF")
    {
      GL2PS_TYPE=GL2PS_PDF;
      GL2PS_EXT=".pdf";
    }
  if (gl2ps_type=="PS")
    {
      GL2PS_TYPE=GL2PS_PS;
      GL2PS_EXT=".ps";
    }
  if (gl2ps_type=="EPS")
    {
      GL2PS_TYPE=GL2PS_EPS;
      GL2PS_EXT=".eps";
    }
  if (gl2ps_type=="SVG")
    {
      GL2PS_TYPE=GL2PS_SVG;
      GL2PS_EXT=".svg";
    }
  if (gl2ps_type=="TEX")
    {
      GL2PS_TYPE=GL2PS_TEX;
      GL2PS_EXT=".tex";
    }

  map<int, unsigned int> DEFAULT_ATTRIB_;
  DEFAULT_ATTRIB_[-1] = GLOBAL_DEFAULT_ATTRIB;
  map <int, vector <string> > SHOW_ONLY_ATOMS_;
  SHOW_ONLY_ATOMS_[-1]= vector <string> (1,"CA");

  unsigned int MS=DEFAULT_ATTRIB_INDEX.size();
  if (MS > DEFAULT_ATTRIB.size())
    MS = DEFAULT_ATTRIB.size();
  for (unsigned int i=0; i < MS;++i)
    DEFAULT_ATTRIB_[DEFAULT_ATTRIB_INDEX[i]] = DEFAULT_ATTRIB[i];


  for (unsigned int i=0; i < SHOW_ATOMS_INDEX.size();++i)
    {
      vector <string> v_tmp;
      for (map <unsigned int, string>::iterator itr=ATOM_NAMES[i].begin();
	   itr != ATOM_NAMES[i].end();
	   ++itr)
	v_tmp.push_back(itr->second);
      if (!v_tmp.empty())
	SHOW_ONLY_ATOMS_[SHOW_ATOMS_INDEX[i]] = v_tmp;
    }



  map <int, vector < pair <unsigned int, vector <unsigned int> > > > ATTRIB_CHAIN_;

  ATTRIB_CHAIN_[-1].push_back(make_pair(GLOBAL_ATTRIB_CHAIN_ATTRIB, GLOBAL_ATTRIB_CHAIN_POS));
  MS = ATTRIB_CHAIN_INDEX.size();
  if (MS > ATTRIB_CHAIN_ATTRIB.size())
    MS = ATTRIB_CHAIN_ATTRIB.size();
  if (MS > ATTRIB_CHAIN_POS.size())
    MS = ATTRIB_CHAIN_POS.size();

  for (unsigned int i=0; i < MS; i++)
    {
      vector <unsigned int> v_tmp;
      for(map<unsigned int, unsigned int>::iterator itr=ATTRIB_CHAIN_POS[i].begin();
	  itr != ATTRIB_CHAIN_POS[i].end();
	  ++itr)
	v_tmp.push_back(itr->second);
      ATTRIB_CHAIN_[ATTRIB_CHAIN_INDEX[i]].push_back(make_pair(ATTRIB_CHAIN_ATTRIB[i], v_tmp  ));
    }

  map <int, bool> VISIBLE_;
  VISIBLE_[-1]=GLOBAL_VISIBLE_;

  MS = VISIBLE_INDEX.size();
  if (MS > VISIBLE_FLAG.size())
    MS = VISIBLE_FLAG.size();

  for (unsigned int i=0; i < MS;++i)
    VISIBLE_[VISIBLE_INDEX[i]] = VISIBLE_FLAG[i];


  if(DEBUG_)
    {
      for(map<int,bool>::iterator itr=VISIBLE_.begin();
	  itr!=VISIBLE_.end();
	  ++itr)
	cout << itr->first << " " << itr->second << endl;
    }
  

  for (unsigned int i=0; i <M.size();++i)
    M[i].set_attrib(DEFAULT_ATTRIB_[-1]);

  for (unsigned int i=0; i <M.size();++i)
    if (DEFAULT_ATTRIB_.find(i) != DEFAULT_ATTRIB_.end())
      M[i].set_attrib(DEFAULT_ATTRIB_[i]);



  for (unsigned int i=0; i <M.size();++i)
    if (ATTRIB_CHAIN_.find(-1) != ATTRIB_CHAIN_.end())
      for(unsigned int j=0; j < ATTRIB_CHAIN_[-1].size();++j )
	M[i].set_attrib_ch(ATTRIB_CHAIN_[-1][j].second, ATTRIB_CHAIN_[-1][j].first);

  for (unsigned int i=0; i <M.size();++i)
    if (ATTRIB_CHAIN_.find(i) != ATTRIB_CHAIN_.end())
      for(unsigned int j=0; j < ATTRIB_CHAIN_[i].size();++j )
	M[i].set_attrib_ch(ATTRIB_CHAIN_[i][j].second, ATTRIB_CHAIN_[i][j].first);


  for (unsigned int i=0; i <M.size();++i)
    if (SHOW_ONLY_ATOMS_[-1][0] != "ALL_ATOMS")
      M[i].show_only(SHOW_ONLY_ATOMS_[-1]);

  for (unsigned int i=0; i <M.size();++i)
    if (SHOW_ONLY_ATOMS_.find(i) != SHOW_ONLY_ATOMS_.end())
      if (SHOW_ONLY_ATOMS_[i][0] != "ALL_ATOMS")
	M[i].show_only(SHOW_ONLY_ATOMS_[i]);


  if (DEBUG_)
    {
      cout << "VISIBLE_[-1] " << VISIBLE_[-1] << endl;
    }
  
  for (unsigned int i=0; i <M.size();++i)
    if ((VISIBLE_.find(-1)!=VISIBLE_.end())
	&& VISIBLE_[-1])
      M[i].show();
    else
      M[i].hide();

  for (unsigned int i=0; i <M.size();++i)
    if (VISIBLE_.find(i) != VISIBLE_.end())
      {
	if (VISIBLE_[i])
	  M[i].show();
	else
	  M[i].hide();
      }

  COLOR= new color_map_(COLOR_MAP,color_ 
			(DEFAULT_COLOR[0],
			 DEFAULT_COLOR[1],
			 DEFAULT_COLOR[2],
			 DEFAULT_COLOR[3]));

  /* for (map <unsigned int, map<unsigned int, double> >::iterator itr=COLOR_MAP.begin();
     itr!=COLOR_MAP.end();
     ++itr)
     {
     cout << itr->first << " " <<  itr->second[0] << " " <<  itr->second[1] << " " << itr->second[2] << " " <<itr->second[3] << " " <<  itr->second[4] << endl;
     }
  */
};
//------------------------------------------------------------
//-------------------------------------------------
void update_show(vector<string> str)
{
  for (unsigned int i=0; i <M.size();++i)
    M[i].show_only(str);
};
//--------------------------------------------------------------
void SelectFromMenu(int idCommand)
{
  switch (idCommand)
    {
    case ATOM_CA:
      update_show(str_CA);
      break;
    case ATOM_CA_CB:
      update_show(str_CA_CB);
      break;
    case ATOM_CA_CB_N:
      update_show(str_CA_CB_N);
      break;
    case ATOM_CA_CB_N_O:
      update_show(str_CA_CB_N_O);
      break;
    case CAPTURE_VIEWPORT:
      CaptureViewPort()  ;
      break;
    case ANIMATE:
      aniOn = !aniOn;
      break;
    case MENU_EXIT:
      exit (0);
      break;
    }
  // Almost any menu selection requires a redraw
  glutPostRedisplay();
}
int BuildPopupMenu (void)
{
  int menu, submenu;
  submenu = glutCreateMenu(SelectFromMenu);
  glutAddMenuEntry("C-alpha",ATOM_CA);
  glutAddMenuEntry("C",ATOM_CA_CB);
  glutAddMenuEntry("N",ATOM_CA_CB_N);
  glutAddMenuEntry("O",ATOM_CA_CB_N_O);


  menu = glutCreateMenu (SelectFromMenu);
  glutAddMenuEntry ("Save Viewport\tS", CAPTURE_VIEWPORT);
  glutAddMenuEntry ("Exit\tQ", MENU_EXIT);
  glutAddMenuEntry ("Animate\tSPACE", ANIMATE);
  glutAddSubMenu("RGB Menu",submenu);
  return menu;
}
//--------------------------------------------------------------

int main(int argc, char *argv[])
{
  configfile="config_iview.cfg";
  string pdbfile="";
  string specstring="";


  options_description infor( "Program information");
  infor.add_options()
    ("help,h", "print help message.")
    ("version,V", "print version number");

  options_description usg( "Usage");
  usg.add_options()
    ("pdbfile,f",value<string>(), "pdbfile []")
    ("configfile,c",value<string>(), "config file  [config_iview.cfg]")
    ("specstring,x",value< string >(), "colorization spec string [] ");
  
  options_description outputopt( "Output options");
  outputopt.add_options()
    ("snapshotfile,o",value<string>(), "snapshotfile file prefix [default: random]")
    ("verbose,v",value< bool >(), "verbosity level [default: 1] ");

  options_description desc( "\n### Yet Another PDB Viewer ### (ishanu chattopadhyay 2017)");
  desc.add(infor).add(usg).add(outputopt);

  positional_options_description p;
  variables_map vm;
  if (argc == 1)
    {
      cout << EMPTY_ARG_MESSAGE << endl;
      return 1;
    }
  try
    {
      store(command_line_parser(argc, argv)
	    .options(desc)
	    .run(), vm);
      notify(vm);
    } 
  catch (std::exception &e)
    {
      cout << endl << e.what() 
	   << endl << desc << endl;
      return 1;
    }
  if (vm.count("help"))
    {
      cout << desc << endl;
      return 1;
    }
  if (vm.count("version"))
    {
      cout << VERSION << endl; 
      return 1;
    }

  if (vm.count("snapshotfile"))
    snapshotfile=vm["snapshotfile"].as<string>();
  if (vm.count("pdbfile"))
    pdbfile=vm["pdbfile"].as<string>();
  if (vm.count("configfile"))
    configfile=vm["configfile"].as<string>();
  if (vm.count("specstring"))
    specstring=vm["specstring"].as<string>();
  if (vm.count("verbose"))
    VERBOSE_=vm["verbose"].as<bool>();

  data_reader(pdbfile, "ATOM", M, specstring);

  glutInit(&argc, argv);
  glutInitDisplayMode (GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
  glutInitWindowSize (winW,winH);
  glutCreateWindow (APP_NAME.c_str());
  glutKeyboardFunc (glutKeyboard);
  getMatrix();  
  glutDisplayFunc (glutDisplay);
  glutReshapeFunc (glutResize);
  glutMotionFunc (glutMotion);
  glutMouseFunc (glutMouse);
  // glutIdleFunc (glutIdle); 
  glInit ();
 

  zprSelectionFunc(glutDisplay);     /* Selection mode draw function */
  zprPickFunc(pick);              /* Pick event client callback   */

  //BuildPopupMenu ();
  //glutAttachMenu (GLUT_MIDDLE_BUTTON);

  glutMainLoop();
  return 1;
}
//-------------------------------------------------


void UTIL_::proc_spec_str(string str,
			  vector<unsigned int>& ATTRIB_CHAIN_INDEX,
			  vector<unsigned int>& ATTRIB_CHAIN_ATTRIB,
			  map<unsigned int,
			  map<unsigned int,
			  unsigned int> >& ATTRIB_CHAIN_POS)
{
  if (DEBUG_)
    cout <<  __PRETTY_FUNCTION__ << endl;

  if (str=="")
    return;

  stringstream ss(str);
  string token;
  ATTRIB_CHAIN_INDEX.clear();
  ATTRIB_CHAIN_ATTRIB.clear();
  ATTRIB_CHAIN_POS.clear();
  unsigned int count=0;
  
  while(getline(ss, token, '#')) 
    {
      unsigned int count1=0;
      stringstream ss1(token);
      string token_index,token_pos,token_att;
      string token_tmp;
      unsigned int count_tmp=0;
      while(getline(ss1,token_tmp,':'))
	{
	  if(count_tmp==0)
	    token_index=token_tmp;
	  if(count_tmp==1)
	    token_pos=token_tmp;
	  if(count_tmp==2)
	    token_att=token_tmp;
	  count_tmp++;	  
	}
      stringstream ss2(token_pos);
      unsigned int token_pos_val;
      map<unsigned int,unsigned int> token_map;
      
      while(ss2>>token_pos_val)
	token_map[count1++]=token_pos_val;

      ATTRIB_CHAIN_INDEX.push_back(atoi(token_index.c_str()));
      ATTRIB_CHAIN_ATTRIB.push_back(atoi(token_att.c_str()));
      ATTRIB_CHAIN_POS[count++]=token_map;
    }
  return;
  
};
