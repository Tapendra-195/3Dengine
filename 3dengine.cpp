#include <iostream>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/core.hpp>
#include <opencv2/videoio.hpp>
#include <opencv2/highgui.hpp>
#include <iostream>
#include <stdio.h>
#include <opencv2/features2d.hpp>
#include <fstream>
#include <ctime>

class Timer
{
public:
    Timer() { clock_gettime(CLOCK_REALTIME, &beg_); }

    double elapsed() {
        clock_gettime(CLOCK_REALTIME, &end_);
        return end_.tv_sec - beg_.tv_sec +
	  (end_.tv_nsec - beg_.tv_nsec)/1000000000.;
    }

    void reset() { clock_gettime(CLOCK_REALTIME, &beg_); }

private:
    timespec beg_, end_;
};





//struct vec3d{
// float x,y,z;
//};
double PI = std::acos(-1);
cv::Matx31d light_dir;

class vertex{
  cv::Matx31d v;
  cv::Matx21d tc;
  vertex(){

  }
  vertex(cv::Matx31d v):v(v){}
  vertex(cv::Matx31d v, cv::Matx21d tc):v(v),tc(tc){

  }
};

class zBuffer{
public:
  double z_buffer[400][400];
  zBuffer(){
    for(int i=0;i<400;i++){
      for(int j=0;j<400;j++){
	z_buffer[i][j]=0;//std::numeric_limits<double>::infinity();
;
      }
    }
  }
  void restart(){
    for(int i=0;i<400;i++){
      for(int j=0;j<400;j++){
	z_buffer[i][j]=0;//std::numeric_limits<double>::infinity();
      }
    }
  }
  bool compare(double z, int y,int x){
    if(z>z_buffer[y][x]){
      //std::cout<<"Prev z = "<<z_buffer[i][j]<<" new = "<<z<<std::endl;
      z_buffer[y][x]=z;
      return true;
    }
    return false;
  }
};

zBuffer zb;

class triangle{
  cv::Matx31d n;//normal
public:
  cv::Matx31d pos;
  cv::Matx31d p[3];
  cv::Vec3b color;
  triangle(){}
  triangle(cv::Matx31d p0,cv::Matx31d p1,cv::Matx31d p2){
    p[0]=p0;
    p[1]=p1;
    p[2]=p2;
    n=find_normal();
    color=cv::Vec3b(0,200,0);
  }

  
  void set_color(cv::Vec3b c){
    //    for(int i=0;i<3;i++){
      color=c;
      // }
  }
  cv::Vec3b get_color(){
    return color;
  }
  void set(cv::Matx31d p0,cv::Matx31d p1,cv::Matx31d p2){
    p[0]=p0;
    p[1]=p1;
    p[2]=p2;
    n=find_normal();
  }
  void set_normal(cv::Matx31d normal){
    n=normal;
  }
  void normalize(){
    for(unsigned i=0;i<3;i++){
      if(p[i](2,0)!=0){
	p[i](0,0) /=p[i](2,0);
	p[i](1,0) /=p[i](2,0);
	p[i](2,0) =1.0/p[i](2,0);
      }
      //      else{
      //p[i](2,0) =-1;
      //}
    }
  }

  double get_z(float xi, float yi){
    float fx=10;
    float fy=10;
    float a=pos(0,0);//p[0](0,0);
    float b=pos(1,0);//p[0](1,0);
    float c=pos(2,0);//p[0](2,0);

    double z=(n(0,0)*a+n(1,0)*b+n(2,0)*c)/((n(0,0)*(xi-200)/fx)+(n(1,0)*(yi-200)/fy)+n(2,0));
    return z;
  }
  void scale(float s){
    for(unsigned i=0;i<3;i++){
      p[i](0,0) =(p[i](0,0)-200)*s+200;
      p[i](1,0) =(p[i](1,0)-200)*s+200;
    }
  }
  cv::Matx31d find_normal(){
    cv::Matx31d a = p[2]-p[0];
    cv::Matx31d b = p[1]-p[0];
    //normal is bXa
    double a1 = a(0,0);
    double a2 = a(1,0);
    double a3 = a(2,0);
    double b1 = b(0,0);
    double b2 = b(1,0);
    double b3 = b(2,0);
    return cv::Matx31d(a2*b3-a3*b2,a3*b1-a1*b3,a1*b2-a2*b1);
  }
  cv::Matx31d get_normal(){
    return n;
  }

  void reevaluate_normal(){
    n=find_normal();
  }

  double area(){
    double area;
    float x1 = p[0](0,0);
    float x2 = p[1](0,0);
    float x3 = p[2](0,0);
    float y1 = p[0](1,0);
    float y2 = p[1](1,0);
    float y3 = p[2](1,0);
    return 0.5*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
  }
  
  void shade(cv::Mat& m, const cv::Matx31d& light_dir){
    cv::Matx31d p_ysorted[3];
    for(unsigned i=0;i<3;i++){
      p_ysorted[i]=p[i];
      p_ysorted[i](1,0)=std::ceil(p_ysorted[i](1,0)-0.5);
    }
    
    for(int i=0;i<3;i++){
      int left=i/2;
      int right=(i+1)%2+1;
      if(p_ysorted[left](1,0)>p_ysorted[right](1,0)){
	cv::Matx31d temp = p_ysorted[right];
	p_ysorted[right]=p_ysorted[left];
	p_ysorted[left]=temp;
      }
    }
    
    //flat top
    if(p_ysorted[0](1,0)==p_ysorted[1](1,0)){
      flat_top(p_ysorted[0],p_ysorted[1],p_ysorted[2],m);
    }
    
    //flat bottom
    else if(p_ysorted[1](1,0)==p_ysorted[2](1,0)){
      flat_bottom(p_ysorted[1],p_ysorted[2],p_ysorted[0],m);
    } 

    else{
      cv::Matx31d intersect;
      intersect(1,0)=p_ysorted[1](1,0);
      intersect(2,0)=p_ysorted[1](2,0);
      float x2=p_ysorted[2](0,0);
      float x1=p_ysorted[0](0,0);
      float y2=p_ysorted[2](1,0);
      float y1=p_ysorted[0](1,0);
      float a = p_ysorted[1](1,0);
      
      intersect(0,0)=((x2-x1)*(a-y1)/(y2-y1))+x1;
      flat_bottom(p_ysorted[1],intersect,p_ysorted[0],m);
      flat_top(p_ysorted[1],intersect,p_ysorted[2],m);
    }
    
  }


  void get_zval(float z1, float z2,float us, float u1, float u2, double& zl ){
    double s=(us-u1)/(u2-u1);
    //zl=(1.0/z1)+s*((1.0/z2)-(1.0/z1));
    zl=z1+s*(z2-z1);
    //    if(zl!=0){
    // zl=1.0/zl;
    //}
  }
  
  void flat_top( cv::Matx31d a, cv::Matx31d b, cv::Matx31d c ,cv::Mat &m){
    if(b(0,0)<a(0,0)){
      cv::Matx31d temp = b;
      b=a;
      a=temp;
    }

    float m1=(c(0,0)-a(0,0))/(c(1,0)-a(1,0));
    float m2=(c(0,0)-b(0,0))/(c(1,0)-b(1,0));
    
    for(unsigned j=a(1,0);j<c(1,0);j++){
      double z1,z2;
      get_zval(a(2,0),c(2,0),j,a(1,0),c(1,0),z1);
      get_zval(b(2,0),c(2,0),j,b(1,0),c(1,0),z2);
      
      float x_low=F(m1,c,j);
      float x_high=F(m2,c,j);
      x_low=ceil(x_low-0.5);
      
      for(int i=x_low;i<x_high;i++){
	//float s=float(i-x_low)/(x_high-x_low);
	double z;
	get_zval(z1,z2,i,x_low,x_high,z);;//(1.0/z1)+s*((1.0/z2)-(1.0/z1));
	//	float z = (1.0/a(2,0))+s*((1.0/b(2,0))-(1.0/a(2,0)));
	//	z=1.0/z;
	
	//std::cout<<"z = "<<z<<std::endl;
	if(zb.compare(z,j,i)){
	  //	  std::cout<<"return true "<<std::endl;
	  float dot = get_normal().dot(light_dir)+0.8;
	  dot=(dot>0)?dot/1.8:0;
	  cv::Vec3b col = get_color();
	  for(unsigned k=0;k<3;k++){
	    col[k] *=dot;
	  }
	  m.at<cv::Vec3b>(j,i)=col;//get_color();//cv::Vec3b(0,100,0);
	}
      }
    }
  }
  
  void flat_bottom(cv::Matx31d a, cv::Matx31d b, cv::Matx31d c,cv::Mat& m ){
    if(b(0,0)<a(0,0)){
      cv::Matx31d temp = b;
      b=a;
      a=temp;
    }
    
    float m1=(a(0,0)-c(0,0))/(a(1,0)-c(1,0));
    float m2=(b(0,0)-c(0,0))/(b(1,0)-c(1,0));
    for(unsigned j=c(1,0);j<a(1,0);j++){
      double z1,z2;
      get_zval(c(2,0),a(2,0),j,c(1,0),a(1,0),z1);
      get_zval(c(2,0),b(2,0),j,c(1,0),b(1,0),z2);

      float x_low=F(m1,c,j);
      float x_high=F(m2,c,j);
      x_low=ceil(x_low-0.5);
      
      for(int i=x_low;i<x_high;i++){  //64,224,208
	double z;
	get_zval(z1,z2,i,x_low,x_high,z);;//(1.0/z1)+s*((1.0/z2)-(1.0/z1));
	//z=1.0/z;
	if(zb.compare(z,j,i)){
	  float dot = get_normal().dot(light_dir)+0.8;
	  dot=(dot>0)?dot/1.8:0;
	  cv::Vec3b col = get_color();
	  for(unsigned k=0;k<3;k++){
	    col[k] *=dot;
	  }
	  m.at<cv::Vec3b>(j,i)=col;//get_color();//cv::Vec3b(0,100,0);
	  
	  // m.at<cv::Vec3b>(j,i)=get_color();//cv::Vec3b(0,100,0);//Vec3b(208,244,64);
	}
      }
    }
  }
  
  double F(double m, cv::Matx31d p,double y){
    
    return m*(y-p(1,0))+p(0,0);
  }
  
};
/*
class 2dtriangle: triangle{
  2dtriangle(){

  }
}
*/
class mesh{
public:
  std::vector<triangle> tris;
  mesh(){
  
  };
  mesh( std::vector<triangle> t):tris(t){}
  /* mesh(){
	   tris=w
  */
  void set(std::vector<triangle> t){
    tris=t;
  }
  void push(triangle t){
    tris.push_back(t);
  }
  
};

class cam{
  float f;//f*no_of_vertical_pixels/cm
  float a;//aspect ratio
  float cx;//width/2
  float cy;//height/2
  cv::Matx33d cam_matrix;
  mesh shape;
  cv::Matx31d dir;
public:
  cv::Mat texture;
  float y; //yaw
  float p; //pitch
  //float dy;
  //float dp;
  cv::Matx31d tr;
  cv::Matx33d R;
  cam(){
    f=10;
    a=1;
    cx=200;
    cy=200;
    cam_matrix= cv::Matx33d(f*a, 0, cx,
			    0,  f,  cy,
			    0,  0,   1 );
    y=0;
    p=0;
    calculate_rotation();
    //dy=0;
    //dp=0;
    tr=cv::Matx31d(0,0,0);
  }
  
  //cv::Matx31d get_translation(){
  //    return
  //}
  void set_mesh(mesh mes){
    shape=mes;
  }
  void calculate_rotation(){
    double thy = y*PI/180.;
    double thx = -p*PI/180.;
    //    y+=dy;
    // p+=dp;
    /*    cv::Matx33d r1=cv::Matx33d(cos(thy),  0, -sin(thy),
			       0,         1, 0,
			       sin(thy), 0, cos(thy));
    
    cv::Matx33d r2=cv::Matx33d(1,  0,        0,
			       0, cos(thx), sin(thx),
			       0, -sin(thx), cos(thx));
    */
    R=cv::Matx33d(cos(thy)         ,0       ,-sin(thy),
		  sin(thx)*sin(thy),cos(thx),sin(thx)*cos(thy),
		  cos(thx)*sin(thy),-sin(thx),cos(thx)*cos(thy)
		  );
    //return r;//r2*r1;
  }

  cv::Matx33d get_rotation(bool changed=false){
    if(changed){
      calculate_rotation();
    }
    return R;
  }
  cv::Matx31d get_dir(cv::Matx31d d){
    dir = (get_rotation().t()*d);
    if(d(1,0)==0){
      dir(1,0)=0;///std::sqrt(dir(0,0)*dir(0,0)+dir(1,0)*dir(1,0)+dir(2,0)*dir(0,0));
    }
    else{
      dir(0,0)=0;
      dir(2,0)=0;
    }
    return dir;
  }

  void set_texture(cv::Mat m){
    texture = m;
  }
  /*
  cv::Matx31d get_right_dir(){
    dir = (get_rotation().t()*cv::Matx31d(1,0,0));
    dir(1,0)=0;///std::sqrt(dir(0,0)*dir(0,0)+dir(1,0)*dir(1,0)+dir(2,0)*dir(0,0));
    return dir;
  }
  */
  void find_intersect(cv::Matx31d p, cv::Matx31d n, cv::Matx31d u, cv::Matx31d v, cv::Matx31d& x){
    double t = n.dot(p-u)/n.dot(v-u);
    x=u+t*(v-u);
  }
  
  int clip(const cv::Matx31d n, cv::Matx31d p, triangle t, triangle& t1, triangle& t2){
    t1.pos=t.pos;
    t2.pos=t.pos;
    t1.set_color(t.get_color());
    t2.set_color(t.get_color());
    //    t2.set_normal(t.get_normal());
    
    //for(
    unsigned in[3];
    int n_in=0;
    for(unsigned i=0;i<3;i++){
      cv::Matx31d u=t.p[i]-p;
      float n_dot_u = n(0,0)*u(0,0)+n(1,0)*u(1,0)+n(2,0)*u(2,0);
      if(n_dot_u>=0){in[n_in++]=i;}
      else{in[2-i+n_in]=i;}
    }

    if(n_in==1){
      //find intersection with line from 3rd index to first
      //find intersection with line from second index to first
      //make a triangle with first point as intersection between first index and third, and second point as first index and second and third as first point.
      cv::Matx31d x0;//lower index one
      find_intersect(p,n,t.p[in[0]],t.p[in[2]],x0);

      cv::Matx31d x1;//higher index one
      find_intersect(p,n,t.p[in[0]],t.p[in[1]],x1);

      if(in[0]==2){
	t1.set(t.p[in[0]],x1,x0);
	t1.set_normal(t.get_normal());
	return n_in;
	  }
      else{
	t1.set(t.p[in[0]],x0,x1);
	t1.set_normal(t.get_normal());
	return n_in;
      }
    }
    else if(n_in==2){
      //find intersection between first index and third index
      //find intersection between second index and third
      //make two triangles:
      //first with first index, second index and intersection of second and third
      //second with 
      cv::Matx31d x0;//lower index one
      find_intersect(p,n,t.p[in[2]],t.p[in[0]],x0);
      
      cv::Matx31d x1;//higher index one
      find_intersect(p,n,t.p[in[2]],t.p[in[1]],x1);
      
      if(in[2]==2){
	t1.set(x1,t.p[in[1]],t.p[in[0]]);
	t2.set(x1,t.p[in[0]],x0);
	//	t1.set(t.p[in[1]],t.p[in[0]],x1);
	//t2.set(t.p[in[0]],x0,x1);
	t1.set_normal(t.get_normal());
	t2.set_normal(t.get_normal());
	return n_in;
      }
      else{
	t1.set(t.p[in[0]],t.p[in[1]],x0);
	t2.set(t.p[in[1]],x1,x0);
	t1.set_normal(t.get_normal());
	t2.set_normal(t.get_normal());
	return n_in;
      }
      
    }
    else if(n_in==3){
      t1.set(t.p[0],t.p[1],t.p[2]);
      t1.set_normal(t.get_normal());
      return n_in;
    }
    else{
      return n_in;
    }
    
  }
  
  cv::Mat get_scene(bool r_changed=false){
    zb.restart();
    double invsqrt2=1.0/sqrt(2);
    light_dir=cv::Matx31d(invsqrt2,-invsqrt2,0);
    cv::Mat img = cv::Mat(400,400,CV_8UC3, cv::Scalar(0,0,0)); 
    mesh m = shape;
    mesh m1;
    cv::Matx33d r = get_rotation(r_changed);

    light_dir=r*(light_dir);
    
    //    for(triangle t: m.tris){
    for(unsigned i=0;i<m.tris.size();i++){
      for(unsigned j=0;j<3;j++){
	m.tris[i].p[j] =r*(m.tris[i].p[j]-tr+cv::Matx31d(0,0,5));
      }
      m.tris[i].reevaluate_normal();
      
      //m.tris[i].p[0](2,0) +=5;
      //m.tris[i].p[1](2,0) +=5;
      //m.tris[i].p[2](2,0) +=5;
      //cv::Matx31d p1 = cam_matrix*t.p[1];
      //cv::Matx31d p2 = cam_matrix*t.p[2];
      //triangle t1(p0,p1,p2);
      
      //  m1.push(t1);
    }

    //find normal after here.
    
    //clipping

    //mesh m2;
    std::vector<triangle>tempt;
    for(triangle t: m.tris){
      triangle t1,t2;
      int no_in=clip(cv::Matx31d(0,0,1),cv::Matx31d(0,0,1),t,t1,t2);
      if(no_in==1||no_in==3){
	tempt.push_back(t1);
      }
      if(no_in==2){
	tempt.push_back(t1);
	tempt.push_back(t2);
      }
     
    }

    m.set(tempt);
    
    for(triangle t: m.tris){
    
      cv::Matx31d p0 = cam_matrix*t.p[0];
      cv::Matx31d p1 = cam_matrix*t.p[1];
      cv::Matx31d p2 = cam_matrix*t.p[2];
      triangle t1(p0,p1,p2);
      t1.pos=t.p[0];
      t1.set_normal(t.get_normal());
      //t1.reevaluate_normal();
      t1.set_color(t.get_color());
      m1.push(t1);
    }
    
    for(unsigned i=0;i<m1.tris.size();i++){
      m1.tris[i].normalize();
      m1.tris[i].scale(20);
    }

    
    //clipping
    //    triangle t1,t2;
    tempt.clear();
    for(triangle t: m1.tris){
      triangle t1,t2;
      int no_in=clip(cv::Matx31d(1,0,0),cv::Matx31d(0,0,0),t,t1,t2);
      if(no_in==1||no_in==3){
	tempt.push_back(t1);
      }
      if(no_in==2){
	tempt.push_back(t1);
	tempt.push_back(t2);
      }
     
    }
    m1.set(tempt);

    tempt.clear();
    for(triangle t: m1.tris){
      triangle t1,t2;
      int no_in=clip(cv::Matx31d(-1,0,0),cv::Matx31d(400,0,0),t,t1,t2);
      if(no_in==1||no_in==3){
	tempt.push_back(t1);
      }
      if(no_in==2){
	tempt.push_back(t1);
	tempt.push_back(t2);
      }
     
    }
    m1.set(tempt);

    tempt.clear();
    for(triangle t: m1.tris){
      triangle t1,t2;
      int no_in=clip(cv::Matx31d(0,1,0),cv::Matx31d(0,0,0),t,t1,t2);
      if(no_in==1||no_in==3){
	tempt.push_back(t1);
      }
      if(no_in==2){
	tempt.push_back(t1);
	tempt.push_back(t2);
      }
     
    }
    m1.set(tempt);

    tempt.clear();
    for(triangle t: m1.tris){
      triangle t1,t2;
      int no_in=clip(cv::Matx31d(0,-1,0),cv::Matx31d(0,400,0),t,t1,t2);
      if(no_in==1||no_in==3){
	tempt.push_back(t1);
      }
      if(no_in==2){
	tempt.push_back(t1);
	tempt.push_back(t2);
      }
     
    }
    m1.set(tempt);
    tempt.clear();
    

    
    for(triangle t: m1.tris){
      if(t.get_normal().dot(t.pos)<0){
	t.shade(img,light_dir);
	//if(true){
	//	line(img, cv::Point2f(t.p[0](0,0),t.p[0](1,0)),cv::Point2f(t.p[1](0,0),t.p[1](1,0)),cv::Scalar(0,255,230),1);
	//line(img, cv::Point2f(t.p[0](0,0),t.p[0](1,0)),cv::Point2f(t.p[2](0,0),t.p[2](1,0)),cv::Scalar(0,255,230),1);
	//line(img, cv::Point2f(t.p[1](0,0),t.p[1](1,0)),cv::Point2f(t.p[2](0,0),t.p[2](1,0)),cv::Scalar(0,255,230),1);
      }
    }
    /*
    cv::Mat A(400, 400, CV_8UC1);
    for(unsigned k=0;k<400;k++){
      for(unsigned l=0;l<400;l++){
	if(zb.z_buffer[k][l]<256){
	  A.at<uchar>(k,l)=255-zb.z_buffer[k][l];
	}
	else{
	  A.at<uchar>(k,l)=0;
	}
      }
    }
    //imshow("z buffer",A);
    */
    
    return img;
  }
};

cam c;
cv::Mat scene;
float prex;
float prey;
bool set=false;
bool changed=false;
mesh m;
//Timer tmr;

void onMouse(int event, int x, int y, int flags, void* userdata){
  bool clicked=false;
  if(set==false && event==0){
    prex=200;
    prey=200;
    set=true;
  }

  if(set){
    float delx = x-prex;
    float dely = y-prey;
    if(delx!=0){
      c.y+=(delx/400.)*800;
      prex=x;
      changed=true;
      
    }
    if(dely!=0){
      c.p+=(dely/400.)*85;
      prey=y;
      changed=true; 
    }
  }

  if(changed){
    scene = c.get_scene(changed);
  }
  //tmr.reset();
  changed=false;
  imshow("Scene", scene);
}


int main(){
  cv::Mat txture=imread("texture.jpg",cv::IMREAD_COLOR);
  std::vector<triangle> t = {
			     //F
			     triangle(cv::Matx31d(-1,-1,-1), cv::Matx31d(1,-1,-1), cv::Matx31d(-1,1,-1)),
			     triangle(cv::Matx31d(-1,1,-1), cv::Matx31d(1,-1,-1), cv::Matx31d(1,1,-1)),
			     //B
			     triangle(cv::Matx31d(1,-1,1), cv::Matx31d(-1,-1,1), cv::Matx31d(1,1,1)),
			     triangle(cv::Matx31d(1,1,1), cv::Matx31d(-1,-1,1), cv::Matx31d(-1,1,1)),
			     //L
			     triangle(cv::Matx31d(-1,-1,1), cv::Matx31d(-1,-1,-1), cv::Matx31d(-1,1,1)),
			     triangle(cv::Matx31d(-1,1,1), cv::Matx31d(-1,-1,-1), cv::Matx31d(-1,1,-1)),
			     //R
			     triangle(cv::Matx31d(1,-1,-1), cv::Matx31d(1,-1,1), cv::Matx31d(1,1,-1)),
			     triangle(cv::Matx31d(1,1,-1), cv::Matx31d(1,-1,1), cv::Matx31d(1,1,1)),
			     //U
			     //triangle(cv::Matx31d(-1,-1,1), cv::Matx31d(1,-1,1), cv::Matx31d(-1,-1,-1)),
			     //triangle(cv::Matx31d(-1,-1,-1), cv::Matx31d(1,-1,1), cv::Matx31d(1,-1,-1)),
			     //D
			     //triangle(cv::Matx31d(-1,1,-1), cv::Matx31d(1,1,-1), cv::Matx31d(-1,1,1)),
			     //triangle(cv::Matx31d(-1,1,1), cv::Matx31d(1,1,-1), cv::Matx31d(1,1,1))
			     
  };
  /*
  for(unsigned i=0;i<t.size();i++){
    t[i].set_texture_coord(cv::Matx21d(txture),cv::Matx21d(),cv::Matx21d());
  }
  */
  
  std::vector<cv::Vec3b> cl ={cv::Vec3b(0,0,200),cv::Vec3b(0,0,200),cv::Vec3b(0,200,0),cv::Vec3b(0,200,0),cv::Vec3b(200,0,0),cv::Vec3b(200,0,0)};
  unsigned col_ind=0;
  for(int i=0;i<t.size();i++){
    cv::Vec3b col =cl[col_ind];
    t[i].set_color(col);
    col_ind++;
    if(col_ind==6){col_ind=0;}
    
  }
  
  int max=t.size();
  for(int i=0;i<max;i++){
    cv::Matx31d ad(4,0,0);
    triangle tr(t[i].p[0]+ad,t[i].p[1]+ad,t[i].p[2]+ad);
    //
    t.push_back(tr);
    ad=cv::Matx31d(4,0,-4);
    tr=triangle(t[i].p[0]-ad,t[i].p[1]-ad,t[i].p[2]-ad);
    tr.set_color(cv::Vec3b(200,0,0));
    t.push_back(tr);
    ad=cv::Matx31d(0,2,0);
    tr=triangle(t[i].p[0]-ad,t[i].p[1]-ad,t[i].p[2]-ad);
    t.push_back(tr);
    ad=cv::Matx31d(0,4,0);
    tr=triangle(t[i].p[0]-ad,t[i].p[1]-ad,t[i].p[2]-ad);
    t.push_back(tr);
    ad=cv::Matx31d(4,2,-4);
    tr=triangle(t[i].p[0]-ad,t[i].p[1]-ad,t[i].p[2]-ad);
    t.push_back(tr);
    /*    ad=cv::Matx31d(4,0,-1);
    tr=triangle(t[i].p[0]-ad,t[i].p[1]-ad,t[i].p[2]-ad);
    t.push_back(tr);*/
  } 

  t.push_back(triangle(cv::Matx31d(-1,-5,1), cv::Matx31d(1,-5,1), cv::Matx31d(-1,-5,-1)));
  t.push_back(triangle(cv::Matx31d(-1,-5,-1), cv::Matx31d(1,-5,1), cv::Matx31d(1,-5,-1)));
  t.push_back(triangle(cv::Matx31d(3,-1,1), cv::Matx31d(5,-1,1), cv::Matx31d(3,-1,-1)));
  t.push_back(triangle(cv::Matx31d(3,-1,-1), cv::Matx31d(5,-1,1), cv::Matx31d(5,-1,-1)));
  t.push_back(triangle(cv::Matx31d(-5,-3,5), cv::Matx31d(-3,-3,5), cv::Matx31d(-5,-3,3)));
  t.push_back(triangle(cv::Matx31d(-5,-3,3), cv::Matx31d(-3,-3,5), cv::Matx31d(-3,-3,3)));
			   

  
  for(int i=-8;i<9;i++){
    for(int j=-8;j<9;j++){
      //if(i==-1&&j==-1)continue;
      triangle tr=triangle(cv::Matx31d(i,1.1,j),cv::Matx31d(i,1.1,j+1),cv::Matx31d(i+1,1.1,j+1));
      triangle tr1=triangle(cv::Matx31d(i,1.1,j),cv::Matx31d(i+1,1.1,j+1),cv::Matx31d(i+1,1.1,j));
      //triangle tr=triangle(cv::Matx31d(i,1,j),cv::Matx31d(i,1,j+1),cv::Matx31d(i+1,1,j+1));
      //triangle tr1=triangle(cv::Matx31d(i,1,j),cv::Matx31d(i+1,1,j+1),cv::Matx31d(i+1,1,j));

      tr.set_color(cv::Vec3b(0,0,200));
      tr1.set_color(cv::Vec3b(0,0,200));
      t.push_back(tr);
      t.push_back(tr1);
    }
  }
  
  m.set(t);
  c.set_mesh(m);
  c.set_texture(txture);
  
 scene = c.get_scene();
 while(true){
   cv::namedWindow("Scene", cv::WINDOW_NORMAL);    
   //cv::resizeWindow("Scene", 800, 800);
   cv::setMouseCallback("Scene", onMouse);
   char key = cv::waitKey(0);
   //std::cout<<"waitkey - "<<key<<std::endl;
   if(key=='a'){ c.tr-=0.2*c.get_dir(cv::Matx31d(1,0,0));  scene = c.get_scene(); }
   else if(key=='s'){c.tr-=0.2*c.get_dir(cv::Matx31d(0,0,1));  scene = c.get_scene();}
   else if(key=='d'){c.tr+=0.2*c.get_dir(cv::Matx31d(1,0,0));  scene = c.get_scene();}
   else if(key=='w'){c.tr+=0.2*c.get_dir(cv::Matx31d(0,0,1));  scene = c.get_scene();}
   
   else if(key==' '){c.tr+=0.2*c.get_dir(cv::Matx31d(0,-1,0));  scene = c.get_scene();}
   else if(key='l'){c.tr+=0.2*c.get_dir(cv::Matx31d(0,1,0));  scene = c.get_scene(); }
   //std::cout<<"tr = "<<c.tr<<std::endl;
   imshow("Scene",scene);
 }  
 
 
 return 0;
}
