package BC;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import Common.Customer;
import Common.Cuts;
import Common.Facility;
import Common.Instance;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloModeler;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.LazyConstraintCallback;

public class OA_gap {
	public static double CalculateGradient(Instance instance,int j,double[] u,double[] v,double[] h,double beta) throws IloException {
		double gradient_j=0;
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faL0=instance.getFL();
		ArrayList<Facility> faF0=instance.getFF();
		ArrayList<Facility> faC=instance.getFcandidate();
		for(int i=0;i<cus.size();i++) {
			double wij=h[i]/calDistanceSquare(cus.get(i),faC.get(j),beta);
			double shareL=0;
			for(int k=0;k<faL0.size();k++) {
				shareL+=faL0.get(k).getSize()/calDistanceSquare(cus.get(i),faL0.get(k),beta);
			}
			for(int k=0;k<faC.size();k++) {
				shareL+=u[k]/calDistanceSquare(cus.get(i),faC.get(k),beta);
			}
			double shareF=0;
			for(int k=0;k<faF0.size();k++) {
				
				shareF+=faF0.get(k).getSize()/calDistanceSquare(cus.get(i),faF0.get(k),beta);
			}
			for(int k=0;k<faC.size();k++) {
				shareF+=v[k]/calDistanceSquare(cus.get(i),faC.get(k),beta);
			}
			gradient_j+=shareL*wij/(shareL+shareF)/(shareL+shareF);
			//System.out.println(gradient_j);
		}
		return gradient_j;
	}
	
	public static double calDistanceSquare(Customer cus,Facility fa,double beta) {
		double ds=0;
		ds=Math.sqrt((cus.getX()-fa.getX())*(cus.getX()-fa.getX())+(cus.getY()-fa.getY())*(cus.getY()-fa.getY()));	
		double ds1=0;
		ds1=Math.pow(ds, beta);
		return ds1;
	}
	
//	public static double calDistanceSquare(Customer cus,Facility fa) {
//		double ds=0;
//		ds=(cus.getX()-fa.getX())*(cus.getX()-fa.getX())+(cus.getY()-fa.getY())*(cus.getY()-fa.getY());	
//		return ds;
//	}
	
	
	public static double calMarketShare(double[] x,double[] y,Instance instance,double[] u,double[] v,double[] h,double beta) {
		double share=0;
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faL0=instance.getFL();
		ArrayList<Facility> faF0=instance.getFF();
		ArrayList<Facility> faC=instance.getFcandidate();
		for(int i=0;i<cus.size();i++) {
			double shareL=0;
			for(int k=0;k<faL0.size();k++) {
				shareL+=faL0.get(k).getSize()/calDistanceSquare(cus.get(i),faL0.get(k),beta);
			}
			for(int k=0;k<faC.size();k++) {
				shareL+=u[k]/calDistanceSquare(cus.get(i),faC.get(k),beta);
			}
			double shareF=0;
			for(int k=0;k<faF0.size();k++) {
				shareF+=faF0.get(k).getSize()/calDistanceSquare(cus.get(i),faF0.get(k),beta);
			}
			for(int k=0;k<faC.size();k++) {
				shareF+=v[k]/calDistanceSquare(cus.get(i),faC.get(k),beta);
			}
			share+=h[i]*shareF/(shareL+shareF);

			
		}
		for(int j=0;j<faC.size();j++) {
			share-=x[j]*y[j];
		}
		return share;
	}
	
	public static Cuts makecuts(IloNumVar[] y, double[][] LSol,double[] ySol, IloModeler ilcplex,IloNumVar theta, Instance instance,double theta0,IloNumVar[] v,double[] vSol,double[] h,double beta) throws IloException {
		Cuts cuts=new Cuts();
		ArrayList<IloNumExpr> cutLhs = new ArrayList<IloNumExpr>();
		ArrayList<Double> cutRhs = new ArrayList<Double>();
		
		IloLinearNumExpr lhsExpr=ilcplex.linearNumExpr();
		double rhsExpr=0;
		double lhs = 0, rhs = 0;
		ArrayList<Facility> faC=instance.getFcandidate();
		
		lhsExpr.addTerm(1, theta); 
		
		rhsExpr+=calMarketShare(LSol[0],ySol,instance,LSol[1],vSol,h,beta);
		
		
		for(int i=0;i<faC.size();i++) {
			lhsExpr.addTerm(-CalculateGradient(instance,i,LSol[1],vSol,h,beta), v[i]);
			lhsExpr.addTerm(LSol[0][i], y[i]);
			rhsExpr-=CalculateGradient(instance,i,LSol[1],vSol,h,beta)*vSol[i];
			rhsExpr+=LSol[0][i]*ySol[i];
		}
		lhs=theta0;
		rhs=calMarketShare(LSol[0],ySol,instance,LSol[1],vSol,h,beta);

		if(lhs>rhs&&(lhs-rhs)>1e-6) {
			cutLhs.add(lhsExpr);
			cutRhs.add(rhsExpr);
		}
		
		cuts.setLhs(cutLhs);
		cuts.setRhs(cutRhs);
		
		return cuts;
	}

	public static double[][] SeparationProblem(double[][] LSol,Instance instance,double V,double[] a,double[] h,int R,double limit,double beta) throws IloException {
		IloCplex ilcplex =new IloCplex();
		ilcplex.setOut(null);
		ilcplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1.0e-1);
		int fnum=instance.getFnum();
	    IloIntVar[] y=new IloIntVar[fnum];  
	    IloNumVar[] v=new IloNumVar[fnum];
	    double[][] sol=new double[2][fnum];
	    for(int i=0;i<fnum;i++) {
	    	y[i]=ilcplex.boolVar();
	    	v[i]=ilcplex.numVar(0, a[i]);
	    }
	    //double LB=calLB(xSol,instance);

	  //objective
	    IloNumVar theta=ilcplex.numVar(0, 1);
	    
	    ilcplex.addMaximize(theta);
	  //constraints
	    
	    IloLinearNumExpr expr1=ilcplex.linearNumExpr();
	    for(int i=0;i<fnum;i++) {
	    	expr1.addTerm(1, y[i]);
	    }
	    ilcplex.addLe(expr1, R);
	    
	    
	    IloLinearNumExpr expr2=ilcplex.linearNumExpr();
	    for(int i=0;i<fnum;i++) {
	    	expr2.addTerm(1, v[i]);
	    }
	    ilcplex.addLe(expr2, V);
	    
	    
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addGe(v[i], ilcplex.prod(limit*a[0], y[i]));
	    }
	    
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addLe(v[i],ilcplex.prod(a[i], y[i]));
	    }
 	    ilcplex.use(new Callback(y,LSol,ilcplex,instance,theta,v,h,beta));
 		ilcplex.use(new LazyCallback(y,LSol,ilcplex,instance,theta,v,h,beta));
 		
 		
	    
	    if(ilcplex.solve()) {
	    	//ilcplex.exportModel("SP.lp");
	    	double[] yVal = ilcplex.getValues(y);
	    	double[] vVal = ilcplex.getValues(v);
	    	sol[0]=yVal;
	    	sol[1]=vVal;
	    }
	    return sol;
	}
	public static class Callback extends IloCplex.UserCutCallback{
		Cuts cut;
		ArrayList<IloNumExpr> cutLhs;
		ArrayList<Double> cutRhs;
		IloIntVar[] y;
		double[][] LSol;
		double[] h;
		IloCplex ilcplex;
		Instance instance;
		IloNumVar theta;
		IloNumVar[] v;
		double beta;
		
		
		Callback(IloIntVar[] y0,double[][] LSol0,IloCplex ilcplex0,Instance instance0,IloNumVar theta0,IloNumVar[] v0,double[] h0,double beta0){
			y=y0;
			LSol=LSol0;
			ilcplex=ilcplex0;
			instance=instance0;
			theta=theta0;
			v=v0;
			h=h0;
			beta=beta0;
		}
		
		public void main() throws IloException{
			//System.out.println("callback");
			double[] ySol =getValues(y);
			double theta0=getValue(theta);
			double[] ySize=getValues(v);
			cut = makecuts(y, LSol,ySol,ilcplex,theta, instance,theta0,v,ySize,h,beta);
		
			cutLhs = cut.getLhs();
			cutRhs= cut.getRhs();
			for(int i = 0; i< cutLhs.size(); i++) {
				addLocal(ilcplex.le(ilcplex.diff(cutLhs.get(i),cutRhs.get(i)), 0));
			}
		}
	}
	
	public static class LazyCallback extends LazyConstraintCallback{
		Cuts cut;
		ArrayList<IloNumExpr> cutLhs;
		ArrayList<Double> cutRhs;
		IloIntVar[] y;
		double[][] LSol;
		IloCplex ilcplex;
		Instance instance;
		IloNumVar theta;
		IloNumVar[] v;
		double[] h;
		double beta;
		
		
		LazyCallback(IloIntVar[] y0,double[][] LSol0,IloCplex ilcplex0,Instance instance0,IloNumVar theta0,IloNumVar[] v0,double[] h0,double beta0){
			y=y0;
			LSol=LSol0;
			ilcplex=ilcplex0;
			instance=instance0;
			theta=theta0;
			v=v0;
			h=h0;
			beta=beta0;
		}
		
		public void main() throws IloException{
			//System.out.println("lazycallback");
			double[] ySol =getValues(y);
			double theta0=getValue(theta);
			double[] ySize=getValues(v);
			cut = makecuts(y, LSol,ySol,ilcplex,theta, instance,theta0,v,ySize,h,beta);
		
			cutLhs = cut.getLhs();
			cutRhs= cut.getRhs();
			for(int i = 0; i< cutLhs.size(); i++) {
				addLocal(ilcplex.le(ilcplex.diff(cutLhs.get(i),cutRhs.get(i)), 0));
			}
		}
		
	}
	public static void main(String[] args) throws IOException, IloException{
		Main m=new Main();
		String pathname2="datasets/Facility_Leader.txt";
		File filename2 = new File(pathname2);
		String pathname3="datasets/Facility_Follower.txt";
		File filename3 = new File(pathname3);
		String pathname1="RandomData/"+"EX1/"+"80-20"+"/"+"80-20"+"-"+1+".txt";
		File filename1 = new File(pathname1);
		Instance instance=new Instance();
		instance=m.initData(filename1,filename2,filename3);
		
		int fnum=instance.getFnum();
		int cnum=instance.getCnum();
		int R=2;
		double V=100;
		double[] h=new double[cnum];
		for(int k=0;k<cnum;k++) {
			h[k]=1.00/cnum;
		}
		double[] a= new double[fnum];
		for(int j=0;j<fnum;j++) {
			a[j]=50;
		}
		//double time1 = System.nanoTime();
		double limit=0.2;
		double beta=2;
		double[][] LeaderSol= {{0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0},{0,0,0,0,0,50,0,0,0,0,0,0,0,0,0,50,0,0,0,0}};
		double[][] FollwerSol=OA_gap.SeparationProblem(LeaderSol,instance,V,a,h,R,limit,beta);
		System.out.println(Arrays.toString(FollwerSol[0]));
		System.out.println(Arrays.toString(FollwerSol[1]));
		
//		double[] x= {0.0, -0.0, 1.0, -0.0, 0.0, 0.0, 0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -0.0, 0.0, 0.0, 0.0, 0.0};
//		double[] y= {0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.0};
//		double[] u= {0.0, 0.0, 30.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//		double[] v= {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 30.0};


//		double m=calMarketShare(LeaderSol[0],FollwerSol[0],instance,LeaderSol[1],FollwerSol[1],h);
//		//double m=calMarketShare(x,y,instance,u,v,h);
//        double time2 = System.nanoTime();                                   
//		double time = (time2 - time1) / 1e9;//求解时间，单位s
//		System.out.println("运行时间"+time);
//		System.out.println(m);
		
//		System.out.println(Arrays.toString(LeaderSol[0]));
//		System.out.println(Arrays.toString(LeaderSol[1]));
//		System.out.println(Arrays.toString(FollwerSol[0]));
//		System.out.println(Arrays.toString(FollwerSol[1]));
		
	}

}
