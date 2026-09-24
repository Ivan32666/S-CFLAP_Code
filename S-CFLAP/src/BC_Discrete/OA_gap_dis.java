package BC_Discrete;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import BC_Discrete.OA_dis.Callback;
import BC_Discrete.OA_dis.LazyCallback;
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

public class OA_gap_dis {
	public static double CalculateGradient(Instance instance,int j,int k,double[][] z,double[][] w,double[] u,double[] h) throws IloException {
		double gradient_j=0;
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faL0=instance.getFL();
		ArrayList<Facility> faF0=instance.getFF();
		ArrayList<Facility> faC=instance.getFcandidate();
		for(int i=0;i<cus.size();i++) {
			double wijk=h[i]*u[k]/calDistanceSquare(cus.get(i),faC.get(j));
			double shareL=0;
			for(int m=0;m<faL0.size();m++) {
				shareL+=faL0.get(m).getSize()/calDistanceSquare(cus.get(i),faL0.get(m));
			}
			for(int m=0;m<faC.size();m++) {
				for(int n=0;n<u.length;n++) {
					shareL+=z[m][n]*u[n]/calDistanceSquare(cus.get(i),faC.get(m));
				}
				
			}
			double shareF=0;
			for(int m=0;m<faF0.size();m++) {
				shareF+=faF0.get(m).getSize()/calDistanceSquare(cus.get(i),faF0.get(m));
			}
			
			for(int m=0;m<faC.size();m++) {
				for(int n=0;n<u.length;n++) {
					shareF+=w[m][n]*u[n]/calDistanceSquare(cus.get(i),faC.get(m));
				}
				
			}
			
			gradient_j+=shareL*wijk/(shareL+shareF)/(shareL+shareF);
			//System.out.println(gradient_j);
		}
		return gradient_j;
	}
	
	
	public static double calDistanceSquare(Customer cus,Facility fa) {
		double ds=0;
		ds=(cus.getX()-fa.getX())*(cus.getX()-fa.getX())+(cus.getY()-fa.getY())*(cus.getY()-fa.getY());	
		return ds;
	}
	
	
	public static double calMarketShare(double[] x,double[] y,Instance instance,double[][] z,double[][] w,double[] u,double[] h) {
		double share=0;
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faL0=instance.getFL();
		ArrayList<Facility> faF0=instance.getFF();
		ArrayList<Facility> faC=instance.getFcandidate();
		for(int i=0;i<cus.size();i++) {
			double shareL=0;
			for(int k=0;k<faL0.size();k++) {
				shareL+=faL0.get(k).getSize()/calDistanceSquare(cus.get(i),faL0.get(k));
			}
			for(int j=0;j<faC.size();j++) {
				for(int k=0;k<u.length;k++) {
					shareL+=z[j][k]*u[k]/calDistanceSquare(cus.get(i),faC.get(j));
				}
				
			}
			double shareF=0;
			for(int k=0;k<faF0.size();k++) {
				shareF+=faF0.get(k).getSize()/calDistanceSquare(cus.get(i),faF0.get(k));
			}
			for(int j=0;j<faC.size();j++) {
				for(int k=0;k<u.length;k++) {
					shareF+=w[j][k]*u[k]/calDistanceSquare(cus.get(i),faC.get(j));
					
				}
			}
			share+=h[i]*shareF/(shareL+shareF);

			
		}
//		for(int j=0;j<faC.size();j++) {
//			share-=x[j]*y[j];
//		}
		return share;
	}
	
	
	public static Cuts makecuts(IloNumVar[] y, Sol LSol,double[] ySol,IloModeler ilcplex,IloNumVar theta, Instance instance,double theta0,IloIntVar[][] w,double[][] wSol,double[] u,double[] h) throws IloException {
		Cuts cuts=new Cuts();
		ArrayList<IloNumExpr> cutLhs = new ArrayList<IloNumExpr>();
		ArrayList<Double> cutRhs = new ArrayList<Double>();
		
		IloLinearNumExpr lhsExpr=ilcplex.linearNumExpr();
		double rhsExpr=0;
		double lhs = 0, rhs = 0;
		ArrayList<Facility> faC=instance.getFcandidate();
		
		lhsExpr.addTerm(1, theta); 
		
		rhsExpr+=calMarketShare(LSol.getLocation(),ySol,instance,LSol.getDesign_discrete(),wSol,u,h);
		
		for(int j=0;j<faC.size();j++) {
			for(int k=0;k<u.length;k++) {
				lhsExpr.addTerm(-CalculateGradient(instance,j,k,LSol.getDesign_discrete(),wSol,u,h), w[j][k]);
				rhsExpr-=CalculateGradient(instance,j,k,LSol.getDesign_discrete(),wSol,u,h)*wSol[j][k];
			}
			
		}
		lhs=theta0;
		rhs=calMarketShare(LSol.getLocation(),ySol,instance,LSol.getDesign_discrete(),wSol,u,h);

		if(lhs>rhs&&(lhs-rhs)>1e-6) {
			cutLhs.add(lhsExpr);
			cutRhs.add(rhsExpr);
		}
		
		cuts.setLhs(cutLhs);
		cuts.setRhs(cutRhs);
		
		return cuts;
	}


	public static Sol SeparationProblem(Sol LSol,Instance instance,double V,double[] u,double[] a,double[] h,int R) throws IloException {
		IloCplex ilcplex =new IloCplex();
		ilcplex.setOut(null);
		ilcplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1.0e-1);
		int fnum=instance.getFnum();
		int K=u.length;
	    IloIntVar[] y=new IloIntVar[fnum];  
	    IloIntVar[][] w=new IloIntVar[fnum][K];
	    Sol sol=new Sol();
	    for(int j=0;j<fnum;j++) {
	    	y[j]=ilcplex.boolVar();
	    	for(int k=0;k<K;k++) {
	    		w[j][k]=ilcplex.boolVar();
	    	}
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
	    	for(int k=0;k<K;k++) {
	    		expr2.addTerm(u[k], w[i][k]);
	    	}
	    }
	    ilcplex.addLe(expr2, V); 
	    
	    for(int j=0;j<fnum;j++) {
	    	IloLinearNumExpr expr4=ilcplex.linearNumExpr();
	    	for(int k=0;k<K;k++) {
		    	expr4.addTerm(1, w[j][k]);
		    }
	    	ilcplex.addEq(expr4,y[j]);
	    }
	    
 	    ilcplex.use(new Callback(y,LSol,ilcplex,instance,theta,w,u,h));
 		ilcplex.use(new LazyCallback(y,LSol,ilcplex,instance,theta,w,u,h));
 		
	    
	    if(ilcplex.solve()) {
	    	//ilcplex.exportModel("SP.lp");
	    	double[] yVal = ilcplex.getValues(y);
	    	double[][] wVal=new double[fnum][K];
	    	for(int j=0;j<fnum;j++) {
	    		for(int k=0;k<K;k++) {
	    			wVal[j][k]=ilcplex.getValue(w[j][k]);
	    		}
	    	}
	    	sol.setLocation(yVal);
	    	sol.setDesign_discrete(wVal);
	    }
	    return sol;
	}
	public static class Callback extends IloCplex.UserCutCallback{
		Cuts cut;
		ArrayList<IloNumExpr> cutLhs;
		ArrayList<Double> cutRhs;
		IloIntVar[] y;
		Sol LSol;
		double[] h;
		IloCplex ilcplex;
		Instance instance;
		IloNumVar theta;
		IloIntVar[][] w;
		double[] u;
		
		
		Callback(IloIntVar[] y0,Sol LSol0,IloCplex ilcplex0,Instance instance0,IloNumVar theta0,IloIntVar[][] w0,double[] u0,double[] h0){
			y=y0;
			LSol=LSol0;
			ilcplex=ilcplex0;
			instance=instance0;
			theta=theta0;
			w=w0;
			u=u0;
			h=h0;
		}
		
		public void main() throws IloException{
			//System.out.println("callback");
			double[] ySol =getValues(y);
			double theta0=getValue(theta);
			double[][] wSol=new double[ySol.length][u.length];
			for(int j=0;j<ySol.length;j++) {
				for(int k=0;k<u.length;k++) {
					wSol[j][k]=getValue(w[j][k]);
				}
			}
			cut = makecuts(y, LSol,ySol,ilcplex,theta, instance,theta0,w,wSol,u,h);
		
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
		Sol LSol;
		IloCplex ilcplex;
		Instance instance;
		IloNumVar theta;
		IloIntVar[][] w;
		double[] u;
		double[] h;
		
		
		LazyCallback(IloIntVar[] y0,Sol LSol0,IloCplex ilcplex0,Instance instance0,IloNumVar theta0,IloIntVar[][] w0,double[] u0,double[] h0){
			y=y0;
			LSol=LSol0;
			ilcplex=ilcplex0;
			instance=instance0;
			theta=theta0;
			w=w0;
			u=u0;
			h=h0;
		}
		
		public void main() throws IloException{
			//System.out.println("lazycallback");
			double[] ySol =getValues(y);
			double theta0=getValue(theta);
			double[][] wSol=new double[ySol.length][u.length];
			for(int j=0;j<ySol.length;j++) {
				for(int k=0;k<u.length;k++) {
					wSol[j][k]=getValue(w[j][k]);
				}
			}
			cut = makecuts(y, LSol,ySol,ilcplex,theta, instance,theta0,w,wSol,u,h);
		
			cutLhs = cut.getLhs();
			cutRhs= cut.getRhs();
			for(int i = 0; i< cutLhs.size(); i++) {
				addLocal(ilcplex.le(ilcplex.diff(cutLhs.get(i),cutRhs.get(i)), 0));
			}
		}
		
	}
}
