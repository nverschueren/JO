#include<GL/glut.h>
#include<cmath>
#include<armadillo>
using namespace std;
using namespace arma;
/// 1/20/22,  This is a DNS for the  Jump Oscillons model
//-----------------------------------------------------------------------------------------------------
//
//                        GLOBAL VARIABLES
//
//-----------------------------------------------------------------------------------------------------


vec par,k;
mat z, k1,k2,k3,k4,zxx,aux; 

int ret,estro,cont,visual,px3,py3;
bool camina,graba,muestra;
double tme,xp,yp;

//-----------------------------------------------------------------------------------------------------
//
//                        DYNAMICAL SYSTEM
//
//-----------------------------------------------------------------------------------------------------
mat rhs(vec par, mat z);
//-----------------------------------------------------------------------------------------------------
//
//                        WINDOWS
//
//-----------------------------------------------------------------------------------------------------
// Instant solution

void Dibuja(void); 
void Inicializa(void);
void reescala(int x, int y);
void muevelo(int x, int y);
double mapeoy(double umin, double umax, double py,double y);
int mapeox(int px,int Nm1,int x);
void Teclado(unsigned char key,int x, int y);
void tiempo(void);
/*
// parameter space
void Dibuja2(void);
void reescala3(int x, int y);
void mouse(int boton,int apretao,int x, int y);
void Inicializa2(void);
double pixel2valor(int pix,int pmax, double umax, double umin);
double patronmarginal(double theta);
*/
int main(int argc, char **argv)
{

  if(argc==3){
    cout << argv[2] << endl;
    cout << argv[1] << endl;
    par.load(argv[2]); z.load(argv[1]);}
  else{
  system("awk '{print $2}' control_model.dat >aux.dat");
  par.load("aux.dat");
  system("rm aux.dat");
  }
  //cout << par << endl;
  camina=true;
  graba=false;
  
  tme=0;
  estro=par(15);
 
  // cont=0;
  // muestra=false;

  z.set_size(par(0),3);
  zxx.set_size(par(0),3);
  k1.set_size(par(0),3);
  k2.set_size(par(0),3);
  k3.set_size(par(0),3);
  k4.set_size(par(0),3);
  aux.set_size(par(0)+2,3); //plus the number of neighbours you want for your finite differences scheme
  if(argc==1){z.zeros();  z.fill(-0.9540); }
 

  k.set_size(par(0));
  
  k.subvec(0,par(0)/2-1)=regspace(0,par(0)/2-1);
  k.subvec(par(0)/2,k.n_rows-1)=regspace(-par(0)/2,-1);
  k=2*datum::pi*k/par(0)/par(1);
  //  k.save("k.dat",raw_ascii);
	   
  
  
  cout << visual << endl;
  
  
  glutInit( &argc, argv);
  glutInitDisplayMode(GLUT_SINGLE| GLUT_RGBA);//desk pc
 // glutInitDisplayMode(GLUT_DOUBLE| GLUT_RGBA);//laptop

  glutInitWindowPosition(0,0);
  glutInitWindowSize(par(11),par(12));
  glutCreateWindow("Instant solution");
  Inicializa();
  glutDisplayFunc(Dibuja);
  glutReshapeFunc(reescala);
  glutMotionFunc(muevelo);
  glutKeyboardFunc(Teclado);  


  /*
  glutInitWindowSize(200,200);
  glutInitWindowPosition(900,300);
  glutCreateWindow("Espacio de Parametros");
  Inicializa2();
  glutDisplayFunc(Dibuja2);
  glutReshapeFunc(reescala3);
  glutMouseFunc(mouse);
  */

  glutIdleFunc(tiempo);

  glutMainLoop();
  

  return 0;
}

mat rhs(vec par, mat z)
{  aux.rows(1,aux.n_rows-2)=z;
  
  //boundary conditions
  if(par(7)==0){aux.row(0)=aux.row(2); aux.row(aux.n_rows-1)=aux.row(aux.n_rows-3);}//neuman bc
  if(par(7)==1){aux.row(0)=aux.row(aux.n_rows-2);aux.row(aux.n_rows-1)=aux.row(1);}//periodic
  //diffusion
  zxx=(aux.rows(2,aux.n_rows-1)+aux.rows(0,aux.n_rows-3)-2*aux.rows(1,aux.n_rows-2))/par(1)/par(1);
  

  
  //  zxx.col(0)=par(5)*zxx.col(0);  zxx.col(1)=par(6)*zxx.col(1);  zxx.col(2)=par(16)*zxx.col(2);
  if(par(7)==2){
  zxx.col(0)=par(5)*real(ifft(-k%k%fft(z.col(0))));
  zxx.col(1)=par(6)*real(ifft(-k%k%fft(z.col(1))));
  zxx.col(2)=par(16)*real(ifft(-k%k%fft(z.col(2))));
  }
  //reaction
  double kk1,kk3,kk4,tau; kk1=par(9); kk3=par(4);tau=par(3);
  
  aux.col(0).rows(1,aux.n_rows-2)=kk1+2*z.col(0)-z.col(0)%z.col(0)%z.col(0)-kk3*z.col(1)-kk4*z.col(2);
  aux.col(1).rows(1,aux.n_rows-2)=(z.col(0)-z.col(1))/tau;
  aux.col(2).rows(1,aux.n_rows-2)=z.col(0)-z.col(2);
  
  return aux.rows(1,aux.n_rows-2)+zxx;
}


void Dibuja()
{
  glClear(GL_COLOR_BUFFER_BIT);
  //  glPointSize(10);
  if(visual==0){glColor3f(0,0,1);}
  if(visual==1){glColor3f(1,0,0);}
  if(visual==2){glColor3f(0,1,0);}
  
  glLineWidth(3);
  glBegin(GL_LINE_STRIP);
  for(double i=0; i<par(0); i++)
    {glVertex2f(i,z(i,visual));} glEnd();

  glutPostRedisplay();
  glutSwapBuffers();

}

void Inicializa(void)
{
  glClearColor(1,1,1,0); //color del background de la ventana (R,G,B)
  glMatrixMode(GL_PROJECTION);//ni idea que hace
  glLoadIdentity();//ni idea
  gluOrtho2D(0,par(0),par(13),par(14)); //hace el mapeo de la ventana en pixeles a [0:10]x[-1.0:20]

}
void reescala(int x, int y)
{  par(11)=x;par(12)=y;glViewport(0, 0, x, y);
  cout << par(11) << " " << par(12)<<endl; Inicializa();
}
void muevelo(int x, int y)
{
  if(x>=0 && x<par(11) && y>0 && y<par(12)){z(mapeox(par(11),par(0)-1,x),visual)=mapeoy(par(13),par(14),par(12),y);}
  else{cout<< "fuera" << endl;}
}

double mapeoy(double umin, double umax, double py,double y)
{return y*(umin-umax)/py+umax;}
int mapeox(int px,int Nm1,int x)
{return x*Nm1/px;}

void Teclado(unsigned char key,int x, int y)
{
  
  switch(key) 
    {
    case 'p':
      camina=!camina; //pause
      if(camina){cout<< "running.."<<endl;}else{cout<<"paused.."<<endl;}
      break;
    case 32: //pause
      camina=!camina;
      if(camina){cout<< "running.."<<endl;}else{cout<<"paused.."<<endl;}
      break;
    case 'r':
	z.zeros();
	cout << "ras" << endl;
	break;
    case 's'://save
      z.save("campo.dat",arma_ascii);
      par.save("parapara.dat",arma_ascii);
      cout << "state saved!" << endl;
      break;
      
    case 'm':
      estro+=50;
      cout << "estro=" << estro << endl;
      break;
    case 'n':
      if(estro>5){
      estro-=5;

      }else{estro=1;}
      cout << "estro=" << estro << endl;
      break;           
    case 'l'://load
      z.load("campo.dat",arma_ascii);
      par.load("parapara.dat");
      cout << "loading state ..." << endl;
      break;
      
    case 'x':
      visual++;
      visual=visual%3;
      if(visual==0)
	{
	  cout << "showing u " << endl;
	}
      else if(visual==1)
	{
	  cout << "showing v " << endl;
	}
      else if(visual==2)
	{
	  cout << "showing w " << endl;}
      else{ cout << "error" << endl;}
      
      break;
      
    case 'a':
      double mm, ma;
      mm=0; ma=0;
      mm=z.col(visual).min(); ma=z.col(visual).max(); 
      par(13)=mm-.25*abs(mm);
      par(14)=ma+.25*abs(ma);
      cout << "autoscaled!"<< endl;
      Inicializa();
      break;
      

    case 't':
      cout << "---RD Model to observe Jump Oscillons (JO)s----" << endl;
      cout <<  "u_t=Du*nabla**2 u +k1+2*u-u**3-k3*v-k4*w"<< endl;
      cout <<  "v_t=Dv*nabla**2 v+(u-v)/tau" << endl;
      cout << "w_t=Dw*nabla**2*w+u-w " << endl;
      cout << "-----Numerical Parameters---------- " << endl;
      cout << "N:" << par(0) << endl;
      cout << "dx:" << par(1) << endl;
      cout << "L=N*dx:"<< par(0)*par(1)<< endl;
      cout << "dt:" << par(2) << endl;
      //  cout << "esto:" << estro<< endl;
      cout << "-----Model Parameters---------- " << endl;
      cout << "tau:" << par(3) << endl;
      cout << "k3:" << par(4) << endl;
      cout << "Du:" << par(5) << endl;
      cout << "Dv:" << par(6)<< endl;
      cout << "k4:" << par(8)<< endl;
      cout << "k1:" << par(9)<< endl;
      cout << "Dw:" << par(16)<< endl;
      cout << "considering ";  if(par(7)==0){cout << "Neumann ";}else{cout << "periodic"; } cout << " boundary conditions" << endl;
      break;
      
    case 'u':
      cout << "updating parameters!" << endl;
      ret=system("awk '{print $2}' control_model.dat >aux.dat");
      par.load("aux.dat");
      ret=system("rm aux.dat");
      break;

    case 'i':
      cout << "<-----------------KEYS -------------->" << endl;
      cout << "p or space bar: pause/play the simulation " <<endl;
      cout << "r: set all the fields to zero " << endl;
      cout << "s: save state"<< endl;
      cout << "l: load state " <<endl;
      cout << "m/n: (m)ore or less estro " << endl;
      cout << "a: autoscale the windows " <<endl;
      cout << "u: upload (re-read) the parameters file" << endl;
      cout << "t: show the parameter values and equation" << endl;
      cout << "i: watch this information "<< endl;
      cout << "k: save/stop saving the xt diagram in the chosen component xtdata.dat at the local dir "<< endl;
      cout << "o: show/hide time and norm of RHS and values of the solution at the left and right boundaries " << endl;

     break;

    case 'k':
      graba=!graba; //save
      if(graba){cout<< "saving xtdata.."<<endl;}else{cout<<"not saving.."<<endl;}
      break;

    case 'o':
      muestra=!muestra;
      cout << "muestra " << muestra << endl;
      cout << "left:" << z.row(0) << endl;
      cout << "right:" << z.row(z.n_rows-1) << endl;
      
      break;
    case 'y':
      k1.zeros();
      k1.col(0)=real(ifft(-k%k%fft(z.col(0))));

      aux.rows(1,aux.n_rows-2)=z;

      if(par(7)==1){aux.row(0)=aux.row(aux.n_rows-2);aux.row(aux.n_rows-1)=aux.row(1);}//periodic
      //diffusion
      zxx=(aux.rows(2,aux.n_rows-1)+aux.rows(0,aux.n_rows-3)-2*aux.rows(1,aux.n_rows-2))/par(1)/par(1);
      k1.col(1)=zxx.col(0);
      k1.save("ffxx.dat",raw_ascii);      
      break;
      
 
    }
  
  glutPostRedisplay();
}





void tiempo(void)
{

  if(camina)
    {
      cont=0;

      
      while(cont<estro)
	{
	k1=rhs(par,z);
	/*	k2=rhs(par,z+0.5*par(2)*k1);
	k3=rhs(par,z+0.5*par(2)*k2);
	k4=rhs(par,z+par(2)*k3);
	z+=par(2)*(k1+2*k2+2*k3+k4)/6.0;  */ //RK4
	z+=par(2)*k1;//euler
	cont++;
	tme+=par(2);
	}cont=0;
      
         if(muestra)
	   {cout <<"STABILITY: RHS= " << norm(rhs(par,z)) << " time=" << tme << "\r" << flush;}
      
	 
      if(graba)
	{FILE *xxt=fopen("xt.dat","a");
	for(int i=0;i<int(par(0));i++)
	  {fprintf(xxt,"%f\t",z(i,visual));}fprintf(xxt,"\n");fclose(xxt);
	}      
    }//if camina
  
  
  glutPostRedisplay();

}


/*
void Dibuja2(void)
{
 
    glClear(GL_COLOR_BUFFER_BIT);  

    glColor3f(1,0,0); //rojo  
    glBegin(GL_LINE_STRIP);for(double i=0; i<par(17); i+=0.03){glVertex2f(i,patronmarginal(i));}glEnd();

    
    glColor3f(0,0,1); //azul  
  glBegin(GL_LINE_STRIP);for(double i=0; i<thetamax; i+=0.03){glVertex2f(i,biestabilidadmarginal1(i));}glEnd();


    glColor3f(0,0,1); //azul  
  glBegin(GL_LINE_STRIP);for(double i=0; i<thetamax; i+=0.03){glVertex2f(i,biestabilidadmarginal2(i));}glEnd();
    
  glColor3f(0,0,0); //negro
  glPointSize(4);
  glBegin(GL_POINTS);
  glVertex2f(par(4),par(5));
  glEnd();

  //glColor3f(0,0,1); //azul
  glPointSize(8);
  glBegin(GL_POINTS);
  glVertex2f(41.0/30.0,sqrt(1021.0)/30.0);
  glEnd();
  


  glutPostRedisplay();
 glutSwapBuffers();

}

void reescala3(int x, int y)
{
  py3=y;
  px3=x;
  glViewport(0, 0, x, y);
  glClear(GL_COLOR_BUFFER_BIT); 
  glutPostRedisplay();
  Inicializa2();

}

void mouse(int boton,int apretao,int x, int y)
{

if(boton==GLUT_LEFT_BUTTON && apretao==GLUT_DOWN)
  {

    xp=pixel2valor(x,px3,par(17),0);
    yp=pixel2valor(py3-y,py3,par(18),0);
    par(4)=xp;
    par(5)=yp;
    cout << "theta=" << par(4) << ", F=" << par(5) << endl;
  }
  glutPostRedisplay();
}
 void Inicializa2(void)// parameter's space
{
  glClear(GL_COLOR_BUFFER_BIT); //asi se arrgla el problema de la ventana transparente
  glClearColor(1,1,1,0); //color del background de la ventana (R,G,B)
  glMatrixMode(GL_PROJECTION);//ni idea que hace
 glLoadIdentity();//ni idea
 gluOrtho2D(0,par(17),0,par(18)); //hace el mapeo en [0,thetamax]x[0,Fmax]
}

 double pixel2valor(int pix,int pmax, double umax, double umin)
{  return (umax-umin)*pix/pmax+umin;}
double patronmarginal(double theta){  return sqrt(2-2*theta+theta*theta);}
*/
