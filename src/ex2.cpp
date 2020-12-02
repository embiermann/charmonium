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

struct CharmState_new {
  int  N;
  unsigned int  J;
  SIGN P;
  SIGN C;
  string name;
  std::vector<double> v; mass;

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




int main (int argc, char * * argv) {
  vector<CharmState> simulated_state={
        {1,0,MINUS, PLUS,"ηc(1S)",2983.4 },
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
  //
 double max=4; // maximum extent of wave function
 double size=400; // number of steps.
 double delta=max/size; // mesh spacing
 double alpha_s = 0.546;
 double b = 0.143;
 double mc = 1.48;
 double u = mc/2;
 double l = 1;


 Variable r; 

 // pick the first two values
 


  // Generate states:
  int count = 0;
  for (unsigned int L=0;L<=2;L++) {


    // basic structure

        // std::cout << solver.eigenvalues()[1]<<L<<";" << std::endl;



    for (unsigned int S=0;S<=1;S++) {
      int JMIN=abs(int(L-S));
      for (unsigned int J=JMIN;J<=L+S;J++) {
        GENFUNCTION V_basic=-4*alpha_s/(3*r) + b*r;
        // hyperfine
        // GENFUNCTION V_hyp = 0*r;
        // // fine structure
        // GENFUNCTION V_fs = 0*r;
        // total potential
        GENFUNCTION V = V_basic;

        //
        MatrixXd H=MatrixXd::Zero(size,size);
        for (int i=0;i<H.rows();i++) {
          int j=i+1;
          H(i,i)+= 1.0/delta/delta/u;
          if (j>0) H(i,i)+= V(j*delta);
          if (i< H.rows()-1) H(i,i+1) -= 1/2.0/delta/delta/u;
          if (i>0) H(i,i-1) -= 1/2.0/delta/delta/u;
          // plus l item
          GENFUNCTION L_item = L*(L+1)/2.0/r/r/u;
          H(i,i)+= L_item(j*delta);
        }
        SelfAdjointEigenSolver<MatrixXd> solver(H);
        // loop through the vector
        for (unsigned)
    }
  }
  
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
  

  
  for (size_t c=0;c<classify.size();c++) {
    //
    // PLOT THE DATA:
    //
    {
      auto end=partition(state.begin(),state.end(), classify[c]);
      for (auto s=state.begin();s!=end;s++) {
	
	
	
	plotFromPDG.push_back(PFPtr(new PlotFunction1D(FixedConstant(s->mass), RealArg::Gt(c+1.0) && RealArg::Lt(c+2.0))));
	
	PlotFunction1D::Properties prop;
	prop.pen.setStyle(Qt::DotLine);
	prop.pen.setWidth(3.0);
	plotFromPDG.back()->setProperties(prop);
	viewEnergy.add(plotFromPDG.back().get());
	
	textFromPDG.push_back(PTPtr(new PlotText(c+1.0, s->mass+100, QString(s->name.c_str()))));
	viewEnergy.add(textFromPDG.back().get());
	
      }
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

