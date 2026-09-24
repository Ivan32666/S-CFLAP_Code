package BC;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import Common.Customer;
import Common.Cuts;
import Common.Facility;
import Common.Instance;
import Common.Solution;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloModeler;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.LazyConstraintCallback;

public class staticCFLP {
	public static double CalculateGradient(Instance instance,int j,double[] u,double[] v,double[] h) throws IloException {
		double gradient_j=0;
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faC=instance.getFcandidate();
		for(int i=0;i<cus.size();i++) {
			double wij=h[i]/calDistanceSquare(cus.get(i),faC.get(j));
			double shareL=0;

			for(int k=0;k<faC.size();k++) {
				shareL+=u[k]/calDistanceSquare(cus.get(i),faC.get(k));
			}
			double shareF=0;

			for(int k=0;k<faC.size();k++) {
				shareF+=v[k]/calDistanceSquare(cus.get(i),faC.get(k));
			}
			gradient_j+=shareL*wij/(shareL+shareF)/(shareL+shareF);
			//System.out.println(gradient_j);
		}
		return gradient_j;
	}
	
	public static double calDistanceSquare(Customer cus,Facility fa) {
		double ds=0;
		ds=(cus.getX()-fa.getX())*(cus.getX()-fa.getX())+(cus.getY()-fa.getY())*(cus.getY()-fa.getY());	
		return ds;
	}
	
	
	public static double calMarketShare(double[] x,double[] y,Instance instance,double[] u,double[] v,double[] h) {
		double share=0;
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faC=instance.getFcandidate();
		for(int i=0;i<cus.size();i++) {
			double shareL=0;
			for(int k=0;k<faC.size();k++) {
				shareL+=u[k]/calDistanceSquare(cus.get(i),faC.get(k));
			}
			
			double shareF=0;

			for(int k=0;k<faC.size();k++) {
				shareF+=v[k]/calDistanceSquare(cus.get(i),faC.get(k));
			}
			share+=h[i]*shareF/(shareL+shareF);

		}

		return share;
	}
	
	
	public static Cuts makecuts(IloNumVar[] y, double[][] LSol,double[] ySol, IloModeler ilcplex,IloNumVar theta, Instance instance,double theta0,IloNumVar[] v,double[] vSol,double[] h) throws IloException {
		Cuts cuts=new Cuts();
		ArrayList<IloNumExpr> cutLhs = new ArrayList<IloNumExpr>();
		ArrayList<Double> cutRhs = new ArrayList<Double>();
		
		IloLinearNumExpr lhsExpr=ilcplex.linearNumExpr();
		double rhsExpr=0;
		double lhs = 0, rhs = 0;
		ArrayList<Facility> faC=instance.getFcandidate();
		
		lhsExpr.addTerm(1, theta); 
		
		rhsExpr+=calMarketShare(LSol[0],ySol,instance,LSol[1],vSol,h);
		
		
		for(int i=0;i<faC.size();i++) {
			lhsExpr.addTerm(-CalculateGradient(instance,i,LSol[1],vSol,h), v[i]);
			rhsExpr-=CalculateGradient(instance,i,LSol[1],vSol,h)*vSol[i];
		}
		lhs=theta0;
		rhs=calMarketShare(LSol[0],ySol,instance,LSol[1],vSol,h);

		if(lhs>rhs&&(lhs-rhs)>1e-6) {
			cutLhs.add(lhsExpr);
			cutRhs.add(rhsExpr);
		}
		
		cuts.setLhs(cutLhs);
		cuts.setRhs(cutRhs);
		
		return cuts;
	}
	public static double[][] LeaderProblem(Instance instance,double U,double[] a,double[] h,int P,double limit) throws IloException{
		IloCplex ilcplex =new IloCplex();
		ilcplex.setOut(null);
		int fnum=instance.getFnum();
		int cnum=instance.getCnum();
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faC=instance.getFcandidate();
	    IloIntVar[] x=new IloIntVar[fnum];  
	    IloNumVar[] u=new IloNumVar[fnum];
	    double[][] sol=new double[2][fnum];
	    for(int i=0;i<fnum;i++) {
	    	x[i]=ilcplex.boolVar();
	    	u[i]=ilcplex.numVar(0, a[i]);
	    }

	  //objective
	    double[][] dist=new double[cnum][fnum];
	    for(int i=0;i<cnum;i++) {
	    	for(int j=0;j<fnum;j++) {
	    		dist[i][j]=1/calDistanceSquare(cus.get(i),faC.get(j));
	    	}
	    }
	    IloNumExpr obj=ilcplex.numExpr();
	    IloNumExpr[] expr1=new IloNumExpr[cnum];
	    for(int i=0;i<cnum;i++) {
	    	expr1[i]=ilcplex.linearNumExpr();
	    	expr1[i]=ilcplex.scalProd(u, dist[i]);
	    }
	    for(int i=0;i<cnum;i++) {
	    	obj=ilcplex.sum(obj, ilcplex.prod(h[i], expr1[i]));
	    }
	    ilcplex.addMaximize(obj);
	    
	  //constraints
	    
	    IloLinearNumExpr expr2=ilcplex.linearNumExpr();
	    for(int i=0;i<fnum;i++) {
	    	expr2.addTerm(1, x[i]);
	    }
	    ilcplex.addLe(expr2, P);
	    
	    
	    IloLinearNumExpr expr3=ilcplex.linearNumExpr();
	    for(int i=0;i<fnum;i++) {
	    	expr3.addTerm(1, u[i]);
	    }
	    ilcplex.addLe(expr3, U);
   
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addLe(u[i],ilcplex.prod(a[i], x[i]));
	    }
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addGe(u[i], ilcplex.prod(limit*a[0], x[i]));
	    }

	    if(ilcplex.solve()) {
	    	double[] xVal = ilcplex.getValues(x);
	    	double[] uVal = ilcplex.getValues(u);
	    	sol[0]=xVal;
	    	sol[1]=uVal;
	    }
	    ilcplex.end();
	    return sol;
	}
	public static double[][] FollowerProblem(double[][] LSol,Instance instance,double V,double[] a,double[] h,int R,double limit) throws IloException {
		IloCplex ilcplex =new IloCplex();
		ilcplex.setOut(null);
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
 	    //ilcplex.use(new Callback(y,LSol,ilcplex,instance,theta,v,h));
 		ilcplex.use(new LazyCallback(y,LSol,ilcplex,instance,theta,v,h));
 		
 		
	    
	    if(ilcplex.solve()) {
	    	double[] yVal = ilcplex.getValues(y);
	    	double[] vVal = ilcplex.getValues(v);
	    	sol[0]=yVal;
	    	sol[1]=vVal;
	    }
	    ilcplex.end();
	    return sol;
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
		
		
		LazyCallback(IloIntVar[] y0,double[][] LSol0,IloCplex ilcplex0,Instance instance0,IloNumVar theta0,IloNumVar[] v0,double[] h0){
			y=y0;
			LSol=LSol0;
			ilcplex=ilcplex0;
			instance=instance0;
			theta=theta0;
			v=v0;
			h=h0;
		}
		
		public void main() throws IloException{
			//System.out.println("lazycallback");
			double[] ySol =getValues(y);
			double theta0=getValue(theta);
			double[] ySize=getValues(v);
			cut = makecuts(y, LSol,ySol,ilcplex,theta, instance,theta0,v,ySize,h);
		
			cutLhs = cut.getLhs();
			cutRhs= cut.getRhs();
			for(int i = 0; i< cutLhs.size(); i++) {
				addLocal(ilcplex.le(ilcplex.diff(cutLhs.get(i),cutRhs.get(i)), 0));
			}
		}
		
	}
	public static Solution solveStaticCFLP(Instance instance,double[] h,double[] a,int P,int R, double U, double V,double limit) throws IloException {
		Solution solution=new Solution();
		double time1 = System.nanoTime();
		double beta=1;
		double[][] LeaderSol= LeaderProblem(instance,U,a,h,P,limit);
		//double[][] FollowerSol=OA.SeparationProblem(LeaderSol, instance, V, a, h, R,limit);
		double[][] FollowerSol=SOCP.SeparationSOCP(instance, LeaderSol, a, h, V, R, beta,limit);
		double time2 = System.nanoTime();                                   
		double time = (time2 - time1) / 1e9;//求解时间，单位s
		double obj=calMarketShare(LeaderSol[0],FollowerSol[0],instance,LeaderSol[1],FollowerSol[1],h);
		solution.setLeaderLoca(LeaderSol[0]);
		solution.setLeaderAttrac(LeaderSol[1]);
		solution.setFollowerLoca(FollowerSol[0]);
		solution.setFollowerAttrac(FollowerSol[1]);
		solution.setTime(time);
		solution.setObj(1-obj);
		return solution;
	}
	
}
