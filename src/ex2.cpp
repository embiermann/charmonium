#include "Eigen/Dense"
#include "QatPlotWidgets/PlotView.h"
#include "QatPlotWidgets/MultipleViewWindow.h"
#include "QatPlotWidgets/MultipleViewWidget.h"
#include "QatPlotWidgets/CustomRangeDivider.h"
#include "QatPlotting/PlotStream.h"
#include "QatPlotting/RealArg.h"
#include "QatPlotting/PlotFunction1D.h"
#include "QatPlotting/PlotProfile.h"
#include "QatPlotting/PlotText.h"
#include "QatGenericFunctions/Variable.h"
#include "QatGenericFunctions/F1D.h"
#include "QatGenericFunctions/FixedConstant.h"
#include "QatGenericFunctions/GaussQuadratureRule.h"
#include "QatGenericFunctions/CubicSplinePolynomial.h"
#include "QatGenericFunctions/CubicSplinePolynomial.h"
#include "QatGenericFunctions/Variable.h"
#include <QApplication>
#include <QMainWindow>
#include <QToolBar>
#include <QTextEdit>
#include <QAction>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <thread>
#include <algorithm>
using namespace Eigen;
using namespace Genfun;
using namespace std;

#include "QatPlotting/PlotPoint.h"
#include <math.h> 
#include "QatGenericFunctions/Exp.h"
#include "QatGenericFunctions/Theta.h"

typedef unique_ptr<PlotFunction1D> PFPtr;
typedef unique_ptr<PlotProfile>    PPPtr;
typedef unique_ptr<PlotText>       PTPtr;

enum SIGN {PLUS, MINUS};

struct CharmState {
  int  N;
  unsigned int  J;
  SIGN P;
  SIGN C;
  string name;
  double mass;
};




class Is {
public:
  Is(int J, SIGN P, SIGN C):J(J),P(P),C(C){}
  bool operator() (const CharmState & state) const {
    return state.J==J && state.P==P && state.C==C;
  }
private:
  const unsigned int J;
  const SIGN P;
  const SIGN C;
};
vector<CharmState> state={{1,0,MINUS, PLUS,"ηc(1S)",2983.4 },
			  {2,0,MINUS, PLUS,"ηc(2S)",3639.2 },
			  {1,0,PLUS, PLUS,"χc0(1P)",3414.75 },
			  {1,1,PLUS, PLUS,"χc1(1P)",3510.66 },
			  {1,2,PLUS, PLUS,"χc2(1P)",3556.20 },
			  {1,1,PLUS, MINUS,"hc(1P)",3525.38 },
			  {2,0,PLUS, PLUS,"χc0(2P)",3860 },
			  {2,1,PLUS, PLUS,"χc1(2P)",3871.69 },
			  {2,2,PLUS, PLUS,"χc2(2P)",3927.20 },
			  {1,1,MINUS,MINUS,"J/ψ(1S)",3096.9},
			  {1,1,MINUS,MINUS,"ψ(2S)",3686.10},
			  {1,1,MINUS,MINUS,"ψ(3770)",3770},
			  {2,1,MINUS,MINUS,"ψ(4160)",4160}
};

vector<Is>         classify={{0,  MINUS,  PLUS },
			     {1,  MINUS,  MINUS},
			     {1,  PLUS,   MINUS},
			     {0,  PLUS,   PLUS },
			     {1,  PLUS,   PLUS },
			     {2,  PLUS,   PLUS }};

// simulated data
struct CharmState_simu {
  unsigned int N;
  unsigned int L;
  double mass;
  string name;
};


vector<CharmState_simu> state_simu;






int main (int argc, char * * argv) {
  //
 double max=20; // maximum extent of wave function
 double size=400; // number of steps.
 double delta=max/size; // mesh spacing
 double alpha_s = 0.546;
 double b = 0.143;
 double mc = 1.48;
 double u = mc/2;
 double sigma = 1.095;
 Exp exp;
 Theta theta;

 Variable r, ss; 
 // pick the first two values
  // Generate states:

  
  QApplication     app(argc,argv);
  
  MultipleViewWindow window;
  QToolBar *toolBar=window.addToolBar("Tools");
  QAction  *quitAction=toolBar->addAction("Quit");
  
  quitAction->setShortcut(QKeySequence("q"));
  
  QObject::connect(quitAction, SIGNAL(triggered()), &app, SLOT(quit()));

  PlotView viewEnergy;
  {
    PRectF rect;
    rect.setXmin(0);
    rect.setXmax(8);
    rect.setYmin(2800);
    rect.setYmax(4300);
    viewEnergy.setRect(rect);
  }
  window.add(&viewEnergy, "Energy levels");

  
  vector<PFPtr> plotFunctions;
  vector<PlotView *> scanView;

  std::vector<std::thread> thread;


  // Hold on to pointers to functions:
  vector<PPPtr> plotSolutions;



  
  // Hold on to pointers to functions:
  vector<PFPtr> plotFromPDG;

  // Hold on to pointers to text:
  vector<PTPtr> textFromPDG;

  PlotProfile prof;
  {
    PlotProfile::Properties prop;
    prop.symbolSize=10;
    prop.pen.setColor("blue");
    prop.brush.setStyle(Qt::SolidPattern);
    prof.setProperties(prop);
  }

  // calculate the H
  int count = 0;
  for (unsigned int L=0;L<=2;L++) {


    // basic structure
    // std::cout << solver.eigenvalues()[1]<<L<<";" << std::endl;
    for (unsigned int S=0;S<=1;S++) {
      int JMIN=abs(int(L-S));
      for (unsigned int J=JMIN;J<=L+S;J++) {
        GENFUNCTION V_basic=-4*alpha_s/(3*r) + b*r;
        // hyperfine
        GENFUNCTION SS = (ss*ss - (3.0/2.0)*(3.0/2.0) -(1.0/2.0)*(3.0/2.0))/2.0;

        GENFUNCTION V_hyp = (32*M_PI*alpha_s/(9*mc*mc))*pow((sigma/sqrt(M_PI)),3)*exp(-sigma*sigma*r*r)*SS(S);
        // fine structure
                // total potential
        GENFUNCTION V_1 = V_basic + V_hyp;
        // // fine structure
        if (J==L+1) {
          GENFUNCTION V_fs = (2.0*alpha_s/(r*r*r) - b/(2.0*r))*(1.0/mc/mc)*(J*(J+1.0)-L*(L+1.0)-S*(S+1.0))/2.0;
          GENFUNCTION T = (4.0*alpha_s/r/r/r)*(-L/(6.0*(2.0*L+3.0)))*(1.0/mc/mc);
          GENFUNCTION V = V_1 + theta(r-40)*(V_fs + T);

          //
          MatrixXd H=MatrixXd::Zero(size,size);
          for (int i=0;i<H.rows();i++) {
  
            int j=i+1;
            double y = j*delta;
            if (J<L+S) {
              y = y + 0.4;
            }
            H(i,i)+= 1.0/delta/delta/u;
            if (j>0) H(i,i)+= V(y);
            if (i< H.rows()-1) H(i,i+1) -= 1/2.0/delta/delta/u;
            if (i>0) H(i,i-1) -= 1/2.0/delta/delta/u;
            // plus l item
            GENFUNCTION L_item = L*(L+1)/2.0/r/r/u;
            H(i,i)+= L_item(y);
          }
        SelfAdjointEigenSolver<MatrixXd> solver(H);
        for (unsigned int N=1; N<=2; N++) {
          // state_simu[count].push_back({L, N, (2*mc + solver.eigenvalues()[0]) * 1000});
          state_simu.push_back(CharmState_simu());
          state_simu[count].L = L;
          state_simu[count].N = N;
          state_simu[count].mass = (2*mc + solver.eigenvalues()[N-1]) * 1000;
          double p = pow(-1, L+1);
          double c = pow(-1, L+S);
          state_simu[count].name = "P=" + to_string(p) + "; C=" + to_string(c);
          if (J== 0 && p==1 && c ==1) {
              prof.addPoint(1.5+3, state_simu[count].mass);
            }
          else if (J== 0 && p==-1 && c ==1) {
              prof.addPoint(1.5, state_simu[count].mass);
            }
          else if (J== 1 && p==-1 && c ==-1) {
              prof.addPoint(1.5+1, state_simu[count].mass);
            }
          else if (J== 1 && p==1 && c ==-1) {
              prof.addPoint(1.5+2, state_simu[count].mass);
            }
          else if (J== 1 && p==1 && c ==1) {
              prof.addPoint(1.5+4, state_simu[count].mass);
            }
          else if (J== 2 && p==1 && c ==1) {
              prof.addPoint(1.5+5, state_simu[count].mass);
            }
          

          // plot
          viewEnergy.add(&prof);
          // viewEnergy.add(new PlotText(J+1-1.0, state_simu[count].mass+100, QString(state_simu[count].name.c_str())));
          std::cout << state_simu[count].mass << std::endl;
          count ++;

        }
        }
        else if (J==L) {
          GENFUNCTION V_fs = (2.0*alpha_s/(r*r*r) - b/(2.0*r))*(1.0/mc/mc)*(J*(J+1.0)-L*(L+1.0)-S*(S+1.0))/2.0;
          GENFUNCTION T = (4.0*alpha_s/r/r/r)/6.0*(1.0/mc/mc);
          GENFUNCTION V = V_1 + theta(r-40)*(V_fs + T);

          //
          MatrixXd H=MatrixXd::Zero(size,size);
        for (int i=0;i<H.rows();i++) {

          int j=i+1;
          double y = j*delta;
          if (J<L+S) {
            y = y + 0.4;
          }
          H(i,i)+= 1.0/delta/delta/u;
          if (j>0) H(i,i)+= V(y);
          if (i< H.rows()-1) H(i,i+1) -= 1/2.0/delta/delta/u;
          if (i>0) H(i,i-1) -= 1/2.0/delta/delta/u;
          // plus l item
          GENFUNCTION L_item = L*(L+1)/2.0/r/r/u;
          H(i,i)+= L_item(j*delta);
        }
        SelfAdjointEigenSolver<MatrixXd> solver(H);
        for (unsigned int N=1; N<=2; N++) {
          // state_simu[count].push_back({L, N, (2*mc + solver.eigenvalues()[0]) * 1000});
          state_simu.push_back(CharmState_simu());
          state_simu[count].L = L;
          state_simu[count].N = N;
          state_simu[count].mass = (2*mc + solver.eigenvalues()[N-1]) * 1000;
          double p = pow(-1, L+1);
          double c = pow(-1, L+S);
          state_simu[count].name = "P=" + to_string(p) + "; C=" + to_string(c);          
          if (J== 0 && p==1 && c ==1) {
              prof.addPoint(1.5+3, state_simu[count].mass);
            }
          else if (J== 0 && p==-1 && c ==1) {
              prof.addPoint(1.5, state_simu[count].mass);
            }
          else if (J== 1 && p==-1 && c ==-1) {
              prof.addPoint(1.5+1, state_simu[count].mass);
            }
          else if (J== 1 && p==1 && c ==-1) {
              prof.addPoint(1.5+2, state_simu[count].mass);
            }
          else if (J== 1 && p==1 && c ==1) {
              prof.addPoint(1.5+4, state_simu[count].mass);
            }
          else if (J== 2 && p==1 && c ==1) {
              prof.addPoint(1.5+5, state_simu[count].mass);
            }

          // plot
          viewEnergy.add(&prof);
          // viewEnergy.add(new PlotText(J+1-1.0, state_simu[count].mass+100, QString(state_simu[count].name.c_str())));
          std::cout << "L=" + to_string(L) + "; N=" + to_string(N) + "; s=" + to_string(S) + "; j=" + to_string(J) + "; " << state_simu[count].mass << std::endl;
          count ++;

        }
        }
        else if (J==L-1) {
          GENFUNCTION V_fs = (2.0*alpha_s/(r*r*r) - b/(2.0*r))*(1.0/mc/mc)*(J*(J+1.0)-L*(L+1.0)-S*(S+1.0))/2.0;
          GENFUNCTION T = (4.0*alpha_s/r/r/r)*((L+1.0)/(6.0*(2.0*L-1.0)))*(1.0/mc/mc);
          GENFUNCTION V = V_1 + theta(r-40)*(V_fs + T);
          //
                  MatrixXd H=MatrixXd::Zero(size,size);
        for (int i=0;i<H.rows();i++) {

          int j=i+1;
          double y = j*delta;
          if (J<L+S) {
            y = y + 0.4;
          }
          H(i,i)+= 1.0/delta/delta/u;
          if (j>0) H(i,i)+= V(y);
          if (i< H.rows()-1) H(i,i+1) -= 1/2.0/delta/delta/u;
          if (i>0) H(i,i-1) -= 1/2.0/delta/delta/u;
          // plus l item
          GENFUNCTION L_item = L*(L+1)/2.0/r/r/u;
          H(i,i)+= L_item(j*delta);
        }
        SelfAdjointEigenSolver<MatrixXd> solver(H);
        for (unsigned int N=1; N<=2; N++) {
          // state_simu[count].push_back({L, N, (2*mc + solver.eigenvalues()[0]) * 1000});
          state_simu.push_back(CharmState_simu());
          state_simu[count].L = L;
          state_simu[count].N = N;
          state_simu[count].mass = (2*mc + solver.eigenvalues()[N-1]) * 1000;
          double p = pow(-1, L+1);
          double c = pow(-1, L+S);
          state_simu[count].name = "P=" + to_string(p) + "; C=" + to_string(c);          
          if (J== 0 && p==1 && c ==1) {
              prof.addPoint(1.5+3, state_simu[count].mass);
            }
          else if (J== 0 && p==-1 && c ==1) {
              prof.addPoint(1.5, state_simu[count].mass);
            }
          else if (J== 1 && p==-1 && c ==-1) {
              prof.addPoint(1.5+1, state_simu[count].mass);
            }
          else if (J== 1 && p==1 && c ==-1) {
              prof.addPoint(1.5+2, state_simu[count].mass);
            }
          else if (J== 1 && p==1 && c ==1) {
              prof.addPoint(1.5+4, state_simu[count].mass);
            }
          else if (J== 2 && p==1 && c ==1) {
              prof.addPoint(1.5+5, state_simu[count].mass);
            }         // plot
          viewEnergy.add(&prof);  
          // viewEnergy.add(new PlotText(J+1-1.0, state_simu[count].mass+100, QString(state_simu[count].name.c_str())));
         // viewEnergy.add(new PlotText(J+1-1.0, state_simu[count].mass+100, QString(state_simu[count].name.c_str())));
          std::cout << state_simu[count].mass << std::endl;
          count ++;

        }
        }
        
        
        

        //


      }
    }
  }

  // plot simmulated data
  for (unsigned int i=1; i<7; i++) {
      }
  

  
  // for (size_t c=0;c<classify.size();c++) {
  for (size_t c=0;c<classify.size();c++) {
    //
    // PLOT THE DATA:
    //
    {
      auto end=partition(state.begin(),state.end(), classify[c]);
      for (auto s=state.begin();s!=end;s++) {

      	// plotFromPDG.push_back(PFPtr(new PlotFunction1D(FixedConstant(s->mass), RealArg::Gt(c+1.0) && RealArg::Lt(c+2.0))));
        
      	
        plotFromPDG.push_back(PFPtr(new PlotFunction1D(FixedConstant(s->mass), RealArg::Gt(c+1.0) && RealArg::Lt(c+2.0))));

      	PlotFunction1D::Properties prop;
      	prop.pen.setStyle(Qt::DotLine);
      	prop.pen.setWidth(3.0);
      	plotFromPDG.back()->setProperties(prop);
      	viewEnergy.add(plotFromPDG.back().get());
      	
      	textFromPDG.push_back(PTPtr(new PlotText(c+1.0, s->mass+100, QString(s->name.c_str()))));
      	viewEnergy.add(textFromPDG.back().get());
	
      }

      // for (auto s=state_simu.begin();s!=end;s++) {

      //   plotFromPDG.push_back(PFPtr(new PlotFunction1D(FixedConstant(s->mass), RealArg::Gt(c+1.0) && RealArg::Lt(c+2.0))));
        
      //   PlotFunction1D::Properties prop;
      //   prop.pen.setStyle(Qt::DotLine);
      //   prop.pen.setWidth(3.0);
      //   plotFromPDG.back()->setProperties(prop);
      //   viewEnergy.add(plotFromPDG.back().get());
        
      //   textFromPDG.push_back(PTPtr(new PlotText(c+1.0, s->mass+100, QString(s->name.c_str()))));
      //   viewEnergy.add(textFromPDG.back().get());
  
      // }
    }
  }



  
  CustomRangeDivider divider;
  QTextEdit edit;
  PlotStream stream(&edit);
  stream
    << "J"
    << PlotStream::Super()
    << "PC"
    << PlotStream::Normal()
    << "="
    << PlotStream::EndP();
  divider.add(0.5, edit.document()->clone());
  string rowLabels[6][2]={{"0", "-+"},
			  {"1", "--"},
			  {"1", "+-"},
			  {"0", "++"},
			  {"1", "++"},
			  {"2", "++"}};
  

  for (int i=0;i<6;i++) {
    stream << PlotStream::Clear() 
	   << rowLabels[i][0]
	   << PlotStream::Super()
	   << rowLabels[i][1]
      	   << PlotStream::Normal()
	   << PlotStream::EndP();
    divider.add(i+1.5, edit.document()->clone());
  }
  viewEnergy.setXRangeDivider(&divider);
  
  {
    PlotStream titleStream(viewEnergy.titleTextEdit());
    titleStream << PlotStream::Clear()
		<< PlotStream::Center() 
		<< PlotStream::Family("Sans Serif") 
		<< PlotStream::Size(16)
		<< "Charmonium system: spectrum"
		<< PlotStream::EndP();
    
    
    PlotStream xLabelStream(viewEnergy.xLabelTextEdit());
    xLabelStream << PlotStream::Clear()
		 << PlotStream::Center()
		 << PlotStream::Family("Sans Serif")
		 << PlotStream::Size(16)
		 << PlotStream::EndP();
    

  }

  for (size_t l=0;l<state.size();l++) {
    MultipleViewWidget *mvw=new MultipleViewWidget();
    std::ostringstream label;
    label << "J=" << state[l].J << "/L=" <<"?" <<"/S="<<"?";
    window.add(mvw, label.str());

    // Plot the wavefunctions
    PRectF rect;
    rect.setXmin(0.00);
    rect.setXmax(0.02);
    rect.setYmin(-0.4);
    rect.setYmax( 0.4);
    PlotView *viewBasis=new PlotView(rect);
    mvw->add(viewBasis, "Wave Functions");

    // Plot the potential:
    PlotView *viewPot=new PlotView();
    mvw->add(viewPot, "Potentials");


  }



  window.show();
  app.exec();
  
  return 1;
}

