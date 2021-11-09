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


class zBuffer{
public:
  double z_buffer[400][400];
  zBuffer(){
    for(int i=0;i<400;i++){
      for(int j=0;j<400;j++){
	z_buffer[i][j]=std::numeric_limits<double>::infinity();
;
      }
    }
  }
  void restart(){
    for(int i=0;i<400;i++){
      for(int j=0;j<400;j++){
	z_buffer[i][j]=std::numeric_limits<double>::infinity();
      }
    }
  }
  bool compare(double z, int y,int x){
    if(z<z_buffer[y][x]){
      //std::cout<<"Prev z = "<<z_buffer[i][j]<<" new = "<<z<<std::endl;
      z_buffer[y][x]=z;
      return true;
    }
    return false;
  }
};

zBuffer zb;

class vertex{
public:
  cv::Matx31d v; //vertex
  cv::Matx21d tc; // texture coordinate
  vertex(){
  }
  vertex(cv::Matx31d v):v(v){}
  vertex(double a, double b, double c){
    v(0,0)=a;
    v(1,0)=b;
    v(2,0)=c;
  }
  vertex(cv::Matx31d v, cv::Matx21d tc):v(v),tc(tc){}

  void normalize(){
    if(v(2,0)!=0){
      v(0,0) /=v(2,0);
      v(1,0) /=v(2,0);
      v(2,0)=v(2,0); 
    }
  }

  vertex operator-(const vertex& v1){
    vertex vr(v-v1.v,tc);
    return vr;
  } 

  vertex operator+(const vertex& v1){
    vertex vr(v+v1.v,tc);
    return vr;
  }

  void set_texture_coord(cv::Matx21d t){
    tc=t;
  }
};

cv::Mat texture=imread("texture.png",cv::IMREAD_COLOR);

class triangle{
  cv::Matx31d n;//normal
public:
  cv::Matx31d pos;
  vertex p[3];
  cv::Vec3b color;
  triangle(){}
  triangle(vertex p0,vertex p1,vertex p2){
    p[0]=p0;
    p[1]=p1;
    p[2]=p2;
    color=cv::Vec3b(0,200,0);
    find_normal();
  }

  
  void set_color(cv::Vec3b c){
    color=c;
  }
  
  cv::Vec3b get_color(){
    return color;
  }
  void set(vertex p0,vertex p1,vertex p2){
    p[0]=p0;
    p[1]=p1;
    p[2]=p2;
  }
  void set_normal(cv::Matx31d normal){
    n=normal;
  }
  void normalize(){
    for(unsigned i=0;i<3;i++){
      p[i].normalize();
    }
  }
  
  void scale(float s){
    for(unsigned i=0;i<3;i++){
      p[i].v(0,0) =(p[i].v(0,0)-200)*s+200;
      p[i].v(1,0) =(p[i].v(1,0)-200)*s+200;
    }
  }
  cv::Matx31d find_normal(){
    cv::Matx31d a = p[2].v-p[0].v;
    cv::Matx31d b = p[1].v-p[0].v;
    //normal is bXa
    double a1 = a(0,0);
    double a2 = a(1,0);
    double a3 = a(2,0);
    double b1 = b(0,0);
    double b2 = b(1,0);
    double b3 = b(2,0);
    cv::Matx31d n = cv::Matx31d(a2*b3-a3*b2,a3*b1-a1*b3,a1*b2-a2*b1);
    float modn = sqrt(n(0,0)*n(0,0)+n(1,0)*n(1,0)+n(2,0)*n(2,0));
    return n/modn;
  }
  cv::Matx31d get_normal(){
    return n;
  }

  void reevaluate_normal(){
    n=find_normal();
  }

  void get_u(cv::Matx21d a,cv::Matx21d b,double z1,double z2, double z,cv::Matx21d& u,double s){
    a=a/z1;
    b=b/z2;
    u=a+s*(b-a);
    u *=z;
  }
  
  void get_zval(float z1, float z2,float us, float u1, float u2, double& zl ){
    double s=(us-u1)/(u2-u1);
    zl=(1.0/z1)+s*((1.0/z2)-(1.0/z1));
    zl=1.0/zl;
  }
  
  void shade(cv::Mat& m, const cv::Matx31d& light_dir){
    vertex p_ysorted[3];
    for(unsigned i=0;i<3;i++){
      p_ysorted[i]=p[i];
      p_ysorted[i].v(1,0)=std::ceil(p_ysorted[i].v(1,0)-0.5);
    }
    
    for(int i=0;i<3;i++){
      int left=i/2;
      int right=(i+1)%2+1;
      if(p_ysorted[left].v(1,0)>p_ysorted[right].v(1,0)){
	vertex temp = p_ysorted[right];
	p_ysorted[right]=p_ysorted[left];
	p_ysorted[left]=temp;
      }
    }
    
    //flat top
    if(p_ysorted[0].v(1,0)==p_ysorted[1].v(1,0)){
      flat_top(p_ysorted[0],p_ysorted[1],p_ysorted[2],m);
    }
    
    //flat bottom
    else if(p_ysorted[1].v(1,0)==p_ysorted[2].v(1,0)){
      flat_bottom(p_ysorted[1],p_ysorted[2],p_ysorted[0],m);
    } 

    else{
      vertex intersect;
      intersect.v(1,0)=p_ysorted[1].v(1,0);
      
      float x2=p_ysorted[2].v(0,0);
      float x1=p_ysorted[0].v(0,0);
      float y2=p_ysorted[2].v(1,0);
      float y1=p_ysorted[0].v(1,0);
      float a = p_ysorted[1].v(1,0);

      double s=(a-y1)/(y2-y1);
      double z;
      double z1=p_ysorted[0].v(2,0);
      double z2=p_ysorted[2].v(2,0);
      get_zval(z1,z2,a,y1,y2,z);
      intersect.v(2,0)=z;//p_ysorted[1].v(2,0);
      cv::Matx21d u;
      get_u(p_ysorted[0].tc,p_ysorted[2].tc,z1,z2,z,u,s);
      intersect.tc=u;
      
      intersect.v(0,0)=((x2-x1)*(a-y1)/(y2-y1))+x1;
      flat_bottom(p_ysorted[1],intersect,p_ysorted[0],m);
      flat_top(p_ysorted[1],intersect,p_ysorted[2],m);
    }
    
  }


  
  void flat_top( vertex a, vertex b, vertex c ,cv::Mat &m){
    if(b.v(0,0)<a.v(0,0)){
      vertex temp = b;
      b=a;
      a=temp;
    }

    /*triangle in form
      a______b
       \    /
	\  /
	 \/
	 c
    */
    float m1=(c.v(0,0)-a.v(0,0))/(c.v(1,0)-a.v(1,0));
    float m2=(c.v(0,0)-b.v(0,0))/(c.v(1,0)-b.v(1,0));
    
    for(unsigned j=a.v(1,0);j<c.v(1,0);j++){
      double z1,z2;
      get_zval(a.v(2,0),c.v(2,0),j,a.v(1,0),c.v(1,0),z1);
      get_zval(b.v(2,0),c.v(2,0),j,b.v(1,0),c.v(1,0),z2);
      
      float x_low=F(m1,c.v,j);
      float x_high=F(m2,c.v,j);
      x_low=ceil(x_low-0.5);

      cv::Matx21d u1,u2;
      double s = (j-a.v(1,0))/(c.v(1,0)-a.v(1,0));
      get_u(a.tc,c.tc,a.v(2,0),c.v(2,0),z1,u1,s);
      get_u(b.tc,c.tc,b.v(2,0),c.v(2,0),z2,u2,s);
      
      for(int i=x_low;i<x_high;i++){
	//float s=float(i-x_low)/(x_high-x_low);
	double z;
	get_zval(z1,z2,i,x_low,x_high,z);;//(1.0/z1)+s*((1.0/z2)-(1.0/z1));
	cv::Matx21d u;
	s=(i-x_low)/(x_high-x_low);
	get_u(u1,u2,z1,z2,z,u,s);
	//	u *=z;
	//std::cout<<"u = "<<u<<std::endl;
	//	float z = (1.0/a(2,0))+s*((1.0/b(2,0))-(1.0/a(2,0)));
	//	z=1.0/z;
	
	//std::cout<<"z = "<<z<<std::endl;
	
	cv::Vec3b col =  texture.at<cv::Vec3b>(u(1,0),u(0,0));//get_color();
	if(col==cv::Vec3b(255,255,255))continue;
	if(zb.compare(z,j,i)){
	  //	  std::cout<<"return true "<<std::endl;
	  float dot = get_normal().dot(light_dir)+3;
	  dot=dot/4;
	  //	  std::cout<<"dot = "<<dot<<std::endl;
	  
	  for(unsigned k=0;k<3;k++){
	    col[k] *=dot;
	  }
	  m.at<cv::Vec3b>(j,i)=col;//get_color();//cv::Vec3b(0,100,0);
	}
      }
    }
  }
  
  void flat_bottom(vertex a, vertex b, vertex c,cv::Mat& m){
    if(b.v(0,0)<a.v(0,0)){
      vertex temp = b;
      b=a;
      a=temp;
    }

    /*
      c
     / \
    /   \
   a-----b
     */
    
    float m1=(a.v(0,0)-c.v(0,0))/(a.v(1,0)-c.v(1,0));
    float m2=(b.v(0,0)-c.v(0,0))/(b.v(1,0)-c.v(1,0));
    for(unsigned j=c.v(1,0);j<a.v(1,0);j++){
      double z1,z2;
      get_zval(c.v(2,0),a.v(2,0),j,c.v(1,0),a.v(1,0),z1);
      get_zval(c.v(2,0),b.v(2,0),j,c.v(1,0),b.v(1,0),z2);

      float x_low=F(m1,c.v,j);
      float x_high=F(m2,c.v,j);
      x_low=ceil(x_low-0.5);

      cv::Matx21d u1,u2;
      double s = (j-c.v(1,0))/(a.v(1,0)-c.v(1,0));
      get_u(c.tc,a.tc,c.v(2,0),a.v(2,0),z1,u1,s);
      get_u(c.tc,b.tc,c.v(2,0),b.v(2,0),z2,u2,s);
      
      for(int i=x_low;i<x_high;i++){  //64,224,208
	double z;
	get_zval(z1,z2,i,x_low,x_high,z);;//(1.0/z1)+s*((1.0/z2)-(1.0/z1));

	cv::Matx21d u;
	s=(i-x_low)/(x_high-x_low);
	get_u(u1,u2,z1,z2,z,u,s);
	
	//z=1.0/z;
	cv::Vec3b col = texture.at<cv::Vec3b>(u(1,0),u(0,0));//get_color();
	if(col==cv::Vec3b(255,255,255))continue;
	if(zb.compare(z,j,i)){
	  float dot = get_normal().dot(light_dir)+3;
	  dot=dot/4;
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

class mesh{
public:
  std::vector<triangle> tris;
  mesh(){
  
  };
  mesh( std::vector<triangle> t):tris(t){}

  void set(std::vector<triangle> t){
    tris=t;
  }
  void push(triangle t){
    tris.push_back(t);
  }
  /*
  bool is_empty(){
    if(tris.size()==0){return true;}
    return false;
  }
  */
};

class cam{
  float f;//f*no_of_vertical_pixels/cm
  float a;//aspect ratio
  float cx;//width/2
  float cy;//height/2
  cv::Matx33d cam_matrix;
  mesh shape;
  cv::Matx31d dir;
  mesh player;
  bool walk;
  unsigned step;
public:
  float y; //yaw
  float p; //pitch
  //float dy;
  //float dp;
  cv::Matx31d tr;
  cv::Matx33d R;
  cam(){
    f=20;
    a=0.8;
    cx=200;
    cy=200;
    cam_matrix= cv::Matx33d(f*a, 0, cx,
			    0,  f,  cy,
			    0,  0,   1 );
    y=0;
    p=0;
    

    tr=cv::Matx31d(0,0.2,-5);
    walk=false;
   
    std::vector <triangle> tlayer={triangle(vertex(cv::Matx31d(-0.1,-0.2,0),cv::Matx21d(286,0)), vertex(cv::Matx31d(0.1,-0.2,0),cv::Matx21d(365,0)), vertex(cv::Matx31d(-0.1,1-0.2,0),cv::Matx21d(286,160))), triangle(vertex(cv::Matx31d(-0.1,1-0.2,0),cv::Matx21d(286,160)), vertex(cv::Matx31d(0.1,-0.2,0),cv::Matx21d(365,0)), vertex(cv::Matx31d(0.1,1-0.2,0),cv::Matx21d(365,160)))};  
    
    mesh plyr;
    plyr.set(tlayer);
    set_player(plyr);
   
    calculate_rotation();  
  }
  

  void set_mesh(mesh mes){
    shape=mes;
  }
  void set_player(mesh mes){
    player = mes;
    std::cout<<"setting mesh"<<std::endl;
  }

  
  void calculate_rotation(){
    double thy = y*PI/180.;
    double thx = -p*PI/180.;
    //    y+=dy;
    // p+=dp;
    cv::Matx33d r1=cv::Matx33d(cos(thy),  0, -sin(thy),
			       0,         1, 0,
			       sin(thy), 0, cos(thy));
    
    cv::Matx33d r2=cv::Matx33d(1,  0,        0,
			       0, cos(thx), sin(thx),
			       0, -sin(thx), cos(thx));
        
    /*
    R=cv::Matx33d(cos(thy)         ,0       ,-sin(thy),
		  sin(thx)*sin(thy),cos(thx),sin(thx)*cos(thy),
		  cos(thx)*sin(thy),-sin(thx),cos(thx)*cos(thy)
		  );
		  //return r;//r2*r1;*/
    R=r2*r1;

    
    triangle t1=triangle(vertex(cv::Matx31d(-0.1,-0.2,0),cv::Matx21d(286,0)), vertex(cv::Matx31d(0.1,-0.2,0),cv::Matx21d(365,0)), vertex(cv::Matx31d(-0.1,1-0.2,0),cv::Matx21d(286,160)));
    triangle t2=triangle(vertex(cv::Matx31d(-0.1,1-0.2,0),cv::Matx21d(286,160)), vertex(cv::Matx31d(0.1,-0.2,0),cv::Matx21d(365,0)), vertex(cv::Matx31d(0.1,1-0.2,0),cv::Matx21d(365,160)));  
    std::vector <triangle> tlayer;
    tlayer.push_back(t1);
    tlayer.push_back(t2);


      for(unsigned i=0;i<tlayer.size();i++){
	for(unsigned j=0;j<3;j++){
	
	tlayer[i].p[j].tc =player.tris[i].p[j].tc;
	}
      }
      
      mesh plyr;
      plyr.set(tlayer);
      set_player(plyr);
      
      for(unsigned i=0;i<player.tris.size();i++){
	for(unsigned j=0;j<3;j++){
	  player.tris[i].p[j].v =r2*player.tris[i].p[j].v+cv::Matx31d(0,0,1.7);
	}
      }
      
  }

  void flip_walk(){
    //keep variable
    //std::cout<<"previous = "<<walk<<std::endl;
    if(walk){
      for(unsigned i=0;i<player.tris.size();i++){
      for(unsigned j=0;j<3;j++){
	player.tris[i].p[j].tc(1,0) -=160;
	//	std::cout<<"failed = "<<player.tris[i].p[j].tc(1,0)<<std::endl;
      }

    }
    }
    else{
      for(unsigned i=0;i<player.tris.size();i++){
	for(unsigned j=0;j<3;j++){
	  player.tris[i].p[j].tc(1,0) +=160;
	}
	//      player.tris[i].reevaluate_normal();
      }
    }
    //  
    walk =!(walk);
    
    //std::cout<<"current = "<<walk<<std::endl;
  }
  
  cv::Matx33d get_rotation(bool changed=false){
    if(changed){
      calculate_rotation();
    }
    return R;
  }
  cv::Matx31d get_dir(float x, float y, float z){
    cv::Matx31d d(x,y,z);
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



  void get_zval(float z1, float z2,float us, float u1, float u2, double& zl ){
    double s=(us-u1)/(u2-u1);
    zl=(1.0/z1)+s*((1.0/z2)-(1.0/z1));
    zl=1.0/zl;
  }

  void get_u(cv::Matx21d a,cv::Matx21d b,double z1,double z2, double z,cv::Matx21d& u,double s){
    a=a/z1;
    b=b/z2;
    u=a+s*(b-a);
    u *=z;
  }

  
  void find_intersect(cv::Matx31d p, cv::Matx31d n, vertex v1, vertex v2, vertex& x){
    cv::Matx31d u=v1.v;
    cv::Matx31d v=v2.v;
    double t = n.dot(p-u)/n.dot(v-u);
    x.v=u+t*(v-u);
       
    if(n(2,0)==1){
      x.tc=v1.tc+t*(v2.tc-v1.tc);
    }
    else{
      x.v(2,0)=(1.0/v1.v(2,0))+t*((1.0/v2.v(2,0))-(1.0/v1.v(2,0)));
      x.v(2,0)=1.0/x.v(2,0);
      
      get_u(v1.tc,v2.tc,v1.v(2,0),v2.v(2,0), x.v(2,0), x.tc,t);
    }
   
  }
  
  int clip(const cv::Matx31d n, cv::Matx31d p, triangle t, triangle& t1, triangle& t2){
    t1.pos=t.pos;
    t2.pos=t.pos;
    t1.set_color(t.get_color());
    t2.set_color(t.get_color());
    t1.set_normal(t.get_normal());
    t2.set_normal(t.get_normal());
    
    unsigned in[3];
    int n_in=0;
    for(unsigned i=0;i<3;i++){
      cv::Matx31d u=t.p[i].v-p;
      float n_dot_u = n(0,0)*u(0,0)+n(1,0)*u(1,0)+n(2,0)*u(2,0);
      if(n_dot_u>=0){in[n_in++]=i;}
      else{in[2-i+n_in]=i;}
    }

    if(n_in==1){
      //find intersection with line from 3rd index to first
      //find intersection with line from second index to first
      //make a triangle with first point as intersection between first index and third, and second point as first index and second and third as first point.
      vertex x0;//lower index one
      find_intersect(p,n,t.p[in[0]],t.p[in[2]],x0);

      vertex x1;//higher index one
      find_intersect(p,n,t.p[in[0]],t.p[in[1]],x1);
      
      if(in[0]==2){
	t1.set(t.p[in[0]],x1,x0);
	return n_in;
	  }
      else{
	t1.set(t.p[in[0]],x0,x1);
	return n_in;
      }
    }
    else if(n_in==2){
      //find intersection between first index and third index
      //find intersection between second index and third
      //make two triangles:
      //first with first index, second index and intersection of second and third
      //second with 
      vertex x0;//lower index one
      find_intersect(p,n,t.p[in[2]],t.p[in[0]],x0);
      
      vertex x1;//higher index one
      find_intersect(p,n,t.p[in[2]],t.p[in[1]],x1);
      
      if(in[2]==2){
	t1.set(x1,t.p[in[1]],t.p[in[0]]);
	t2.set(x1,t.p[in[0]],x0);
	return n_in;
      }
      else{
	t1.set(t.p[in[0]],t.p[in[1]],x0);
	t2.set(t.p[in[1]],x1,x0);
	return n_in;
      }
      
    }
    else if(n_in==3){
      t1.set(t.p[0],t.p[1],t.p[2]);
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
    
    for(unsigned i=0;i<m.tris.size();i++){
      for(unsigned j=0;j<3;j++){
	m.tris[i].p[j].v =r*(m.tris[i].p[j].v-tr)+cv::Matx31d(0,0,1.7);
      }
      m.tris[i].reevaluate_normal();
      
    }
    
    
    
    //clipping

    //mesh m2;
    std::vector<triangle>tempt;
    for(triangle t: m.tris){
      triangle t1,t2;
      int no_in=clip(cv::Matx31d(0,0,1),cv::Matx31d(0,0,0.3),t,t1,t2);
      if(no_in==1||no_in==3){
	tempt.push_back(t1);
      }
      else if(no_in==2){
	tempt.push_back(t1);
	tempt.push_back(t2);
      }
     
    }

    m.set(tempt);
    for(unsigned i=0;i<player.tris.size();i++){
      player.tris[i].reevaluate_normal();
      m.tris.push_back(player.tris[i]);
    }
    
    for(triangle t: m.tris){
      cv::Matx31d p0 = cam_matrix*t.p[0].v;
      cv::Matx31d p1 = cam_matrix*t.p[1].v;
      cv::Matx31d p2 = cam_matrix*t.p[2].v;
      
      triangle t1(vertex(p0,t.p[0].tc),vertex(p1,t.p[1].tc),vertex(p2,t.p[2].tc));
      t1.pos=t.p[0].v;
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
      else if(no_in==2){
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
      else if(no_in==2){
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
      else if(no_in==2){
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
      else if(no_in==2){
	tempt.push_back(t1);
	tempt.push_back(t2);
      }
     
    }
    m1.set(tempt);
    tempt.clear();
    

    
    for(triangle t: m1.tris){
      if(t.get_normal().dot(t.pos)<0){
	t.shade(img,light_dir);
	/*	if(false){
	  line(img, cv::Point2f(t.p[0].v(0,0),t.p[0].v(1,0)),cv::Point2f(t.p[1].v(0,0),t.p[1].v(1,0)),cv::Scalar(0,255,230),1);
	  line(img, cv::Point2f(t.p[0].v(0,0),t.p[0].v(1,0)),cv::Point2f(t.p[2].v(0,0),t.p[2].v(1,0)),cv::Scalar(0,255,230),1);
	  line(img, cv::Point2f(t.p[1].v(0,0),t.p[1].v(1,0)),cv::Point2f(t.p[2].v(0,0),t.p[2].v(1,0)),cv::Scalar(0,255,230),1);
	  }*/
      }
    }
    
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
 
  float delx = x-200;//prex;
  float dely = y-200;//prey;
    if(delx!=0){
      c.y=(delx/200.)*400;
      prex=x;
      changed=true;
      
    }
    if(dely!=0){
      c.p=(dely/200.)*28;
      prey=y;
      changed=true; 
    }
    
    if(changed){
      scene = c.get_scene(changed);
      imshow("Scene", scene);
      changed=false;
    }
    
}


int main(){
  std::cout<<"setting up t"<<std::endl;
  prex=200;
  prey=200;
 
  //cv::Mat txture=
  std::vector<triangle> t = {
			     //F
	 		     triangle(vertex(cv::Matx31d(-1,-1,-1),cv::Matx21d(0,100)), vertex(cv::Matx31d(1,-1,-1),cv::Matx21d(261,100)), vertex(cv::Matx31d(-1,1,-1),cv::Matx21d(0,0))),
			     //		     triangle(vertex(-1,1,-1), vertex(1,-1,-1), vertex(1,1,-1)),
			     triangle(vertex(cv::Matx31d(-1,1,-1),cv::Matx21d(0,0)), vertex(cv::Matx31d(1,-1,-1),cv::Matx21d(261,100)), vertex(cv::Matx31d(1,1,-1),cv::Matx21d(261,0))),
			     //B
			     triangle(vertex(cv::Matx31d(1,-1,1),cv::Matx21d(0,100)), vertex(cv::Matx31d(-1,-1,1),cv::Matx21d(261,100)), vertex(cv::Matx31d(1,1,1),cv::Matx21d(0,0))),
			     triangle(vertex(cv::Matx31d(1,1,1),cv::Matx21d(0,0)), vertex(cv::Matx31d(-1,-1,1),cv::Matx21d(261,100)), vertex(cv::Matx31d(-1,1,1),cv::Matx21d(261,0))),
			     //L
			     triangle(vertex(cv::Matx31d(-1,-1,1),cv::Matx21d(0,100)), vertex(cv::Matx31d(-1,-1,-1),cv::Matx21d(261,100)), vertex(cv::Matx31d(-1,1,1),cv::Matx21d(0,0))),
			     triangle(vertex(cv::Matx31d(-1,1,1),cv::Matx21d(0,0)), vertex(cv::Matx31d(-1,-1,-1),cv::Matx21d(261,100)), vertex(cv::Matx31d(-1,1,-1),cv::Matx21d(261,0))),
			     //R
			     triangle(vertex(cv::Matx31d(1,-1,-1),cv::Matx21d(0,100)), vertex(cv::Matx31d(1,-1,1),cv::Matx21d(261,100)), vertex(cv::Matx31d(1,1,-1),cv::Matx21d(0,0))),
			     triangle(vertex(cv::Matx31d(1,1,-1),cv::Matx21d(0,0)), vertex(cv::Matx31d(1,-1,1),cv::Matx21d(261,100)), vertex(cv::Matx31d(1,1,1),cv::Matx21d(261,0))),
			     //U
			     //triangle(vertex(-1,-1,1), vertex(1,-1,1), vertex(-1,-1,-1)),
			     //triangle(vertex(-1,-1,-1), vertex(1,-1,1), vertex(1,-1,-1)),
			     //D
			     //triangle(vertex(-1,1,-1), vertex(1,1,-1), vertex(-1,1,1)),
			     //triangle(vertex(-1,1,1), vertex(1,1,-1), vertex(1,1,1))
			     
			       };
 
  
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
    vertex ad(4,0,0);
    triangle tr(t[i].p[0]+ad,t[i].p[1]+ad,t[i].p[2]+ad);
    //
    t.push_back(tr);
    ad=cv::Matx31d(4,0,-4);
    tr=triangle(t[i].p[0]-ad,t[i].p[1]-ad,t[i].p[2]-ad);
    tr.set_color(cv::Vec3b(200,0,0));
    t.push_back(tr);
    ad=vertex(0,2,0);
    tr=triangle(t[i].p[0]-ad,t[i].p[1]-ad,t[i].p[2]-ad);
    t.push_back(tr);
    ad=vertex(0,4,0);
    tr=triangle(t[i].p[0]-ad,t[i].p[1]-ad,t[i].p[2]-ad);
    t.push_back(tr);
    ad=vertex(4,2,-4);
    tr=triangle(t[i].p[0]-ad,t[i].p[1]-ad,t[i].p[2]-ad);
    t.push_back(tr);
    
    }
  		     
  t.push_back(triangle(vertex(cv::Matx31d(-1,-5,1), cv::Matx21d(0,100)), vertex(cv::Matx31d(1,-5,1),cv::Matx21d(100,100)), vertex(cv::Matx31d(-1,-5,-1),cv::Matx21d(0,0))));
  t.push_back(triangle(vertex(cv::Matx31d(-1,-5,-1),cv::Matx21d(0,0)), vertex(cv::Matx31d(1,-5,1),cv::Matx21d(100,100)), vertex(cv::Matx31d(1,-5,-1),cv::Matx21d(100,0))));
    t.push_back(triangle(vertex(cv::Matx31d(3,-1,1),cv::Matx21d(0,100)), vertex(cv::Matx31d(5,-1,1),cv::Matx21d(100,100)), vertex(cv::Matx31d(3,-1,-1),cv::Matx21d(0,0))));
  t.push_back(triangle(vertex(cv::Matx31d(3,-1,-1),cv::Matx21d(0,0)), vertex(cv::Matx31d(5,-1,1),cv::Matx21d(100,100)), vertex(cv::Matx31d(5,-1,-1),cv::Matx21d(100,0))));
  t.push_back(triangle(vertex(cv::Matx31d(-5,-3,5),cv::Matx21d(0,100)), vertex(cv::Matx31d(-3,-3,5),cv::Matx21d(100,100)), vertex(cv::Matx31d(-5,-3,3),cv::Matx21d(0,0))));
  t.push_back(triangle(vertex(cv::Matx31d(-5,-3,3),cv::Matx21d(0,0)), vertex(cv::Matx31d(-3,-3,5),cv::Matx21d(100,100)), vertex(cv::Matx31d(-3,-3,3),cv::Matx21d(100,0))));
  			   
			     
  
  for(int i=-15;i<15;i++){
    for(int j=-15;j<15;j++){
      //if(i==-1&&j==-1)continue;
      triangle tr=triangle(vertex(cv::Matx31d(i,1,j),cv::Matx21d(0,430)),vertex(cv::Matx31d(i,1,j+1),cv::Matx21d(100,430)),vertex(cv::Matx31d(i+1,1,j+1),cv::Matx21d(0,330)));
      triangle tr1=triangle(vertex(cv::Matx31d(i,1,j),cv::Matx21d(0,330)),vertex(cv::Matx31d(i+1,1,j+1),cv::Matx21d(100,430)),vertex(cv::Matx31d(i+1,1,j),cv::Matx21d(100,330)));
      //triangle tr=triangle(vertex(cv::Matx31d(i,1.1,j),cv::Matx21d(0,300)),vertex(cv::Matx31d(i,1.1,j+1),cv::Matx21d(100,300)),vertex(cv::Matx31d(i+1,1.1,j+1),cv::Matx21d(0,200)));
      //triangle tr1=triangle(vertex(cv::Matx31d(i,1.1,j),cv::Matx21d(0,200)),vertex(cv::Matx31d(i+1,1.1,j+1),cv::Matx21d(100,300)),vertex(cv::Matx31d(i+1,1.1,j),cv::Matx21d(100,200)));
      //triangle tr=triangle(vertex(i,1,j),vertex(i,1,j+1),vertex(i+1,1,j+1));
      //triangle tr1=triangle(vertex(i,1,j),vertex(i+1,1,j+1),vertex(i+1,1,j));

      tr.set_color(cv::Vec3b(0,0,200));
      tr1.set_color(cv::Vec3b(0,0,200));
      t.push_back(tr);
      t.push_back(tr1);
    }
  }
  
  m.set(t);
  c.set_mesh(m);
  //c.set_texture(txture);

  std::cout<<"got after setting mesh = "<<std::endl;
  
 scene = c.get_scene();
 while(true){
   cv::namedWindow("Scene", cv::WINDOW_NORMAL);    
   cv::resizeWindow("Scene", 1200, 800);
   cv::setMouseCallback("Scene", onMouse);
   char key = cv::waitKey(0);
   // std::cout<<"key = "<<key<<std::endl;
   
   if(key=='a'){ c.tr-=0.2*c.get_dir(1,0,0); c.flip_walk(); scene = c.get_scene(); imshow("Scene",scene);}
   else if(key=='s'){c.tr-=0.2*c.get_dir(0,0,1); c.flip_walk(); scene = c.get_scene(); imshow("Scene",scene);}
   else if(key=='d'){c.tr+=0.2*c.get_dir(1,0,0); c.flip_walk(); scene = c.get_scene(); imshow("Scene",scene);}
   else if(key=='w'){c.tr+=0.2*c.get_dir(0,0,1); c.flip_walk(); scene = c.get_scene(); imshow("Scene",scene);}
   
   else if(key==' '){c.tr+=0.2*cv::Matx31d(0,-1,0); scene = c.get_scene(); imshow("Scene",scene);}
   else if(key=='c'){c.tr+=0.2*cv::Matx31d(0,1,0); scene = c.get_scene(); imshow("Scene",scene);}    
   
 }  
 
 
 return 0;
}
