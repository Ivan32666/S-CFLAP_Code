package BC_Discrete;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import Common.Customer;
import Common.Facility;
import Common.Instance;
import ilog.concert.*;
import ilog.cplex.*;


public class SOCP_dis {
	
	public static double calDistanceSquare(Customer cus,Facility fa) {
		double ds=0;
		ds=(cus.getX()-fa.getX())*(cus.getX()-fa.getX())+(cus.getY()-fa.getY())*(cus.getY()-fa.getY());	
		return ds;
	}
//	public static double calDistanceSquare(Customer cus,Facility fa) {
//		double ds=0;
//		ds=Math.sqrt((cus.getX()-fa.getX())*(cus.getX()-fa.getX())+(cus.getY()-fa.getY())*(cus.getY()-fa.getY()));	
//		double ds1=0;
//		ds1=Math.pow(ds, 2);
//		return ds1;
//	}
	
	public static double[] calt(Instance instance,double[][] z ,double[] u){
		int cnum=instance.getCnum();
		double[] t=new double[cnum];
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faC=instance.getFcandidate();
		int K=u.length;
		for(int i=0;i<cnum;i++) {
			for(int j=0;j<faC.size();j++) {
				for(int k=0;k<K;k++) {
					t[i]+=z[j][k]*u[k]/calDistanceSquare(cus.get(i),faC.get(j));
				}
				
			}
		}
		
		return t;
	}
	
	

	public static Sol SeparationSOCP(Instance instance,Sol LeaderSol,double[] u,double[] a,double[] h,double V,int R) throws IloException {
		IloCplex ilcplex =new IloCplex();
		ilcplex.setOut(null);
		int fnum=instance.getFnum();
		int cnum=instance.getCnum();
		int K=u.length;
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faC=instance.getFcandidate();
		
	    IloIntVar[] y=new IloIntVar[fnum];  
	    IloIntVar[][] w=new IloIntVar[fnum][K];
	    IloNumVar[] zeta=new IloNumVar[cnum];  
	    IloNumVar[] b=new IloNumVar[cnum];
	    //IloIntVar[] tl=new IloIntVar[fnum];
	    
	    Sol sol=new Sol();
	    for(int i=0;i<fnum;i++) {
	    	y[i]=ilcplex.boolVar("y"+i);
	    	//tl[i]=ilcplex.boolVar();
	    }
	    for(int j=0;j<fnum;j++) {
	    	for(int k=0;k<K;k++) {
	    		w[j][k]=ilcplex.boolVar("w"+j+k);
	    	}
	    }
	    for(int i=0;i<cnum;i++) {
	    	zeta[i]=ilcplex.numVar(0, 1,"zeta"+i);
	    	b[i]=ilcplex.numVar(1, Double.MAX_VALUE,"b"+i);
	    }

	  //objective
	    IloNumVar omega=ilcplex.numVar(0, 1);
	    
	    ilcplex.addMinimize(omega);
	  //constraints
	    
	    
	    IloLinearNumExpr expr1=ilcplex.linearNumExpr();
	    for(int j=0;j<fnum;j++) {
	    	for(int k=0;k<K;k++) {
	    		expr1.addTerm(u[k], w[j][k]);
	    	}
	    	
	    }
	    ilcplex.addLe(expr1, V);
	    
	    IloLinearNumExpr expr2=ilcplex.linearNumExpr();
	    for(int i=0;i<fnum;i++) {
	    	expr2.addTerm(1, y[i]);
	    }
	    ilcplex.addLe(expr2, R);
	   
	    for(int j=0;j<fnum;j++) {
	    	IloLinearNumExpr expr3=ilcplex.linearNumExpr();
	    	for(int k=0;k<K;k++) {
		    	expr3.addTerm(1, w[j][k]);
		    }
	    	ilcplex.addEq(expr3,y[j]);
	    }
	   
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addLe(y[i],1-LeaderSol.getLocation()[i]);
	    }
	    
	    double[] t=calt(instance,LeaderSol.getDesign_discrete(),u);
	    
	    
	    double[][] dist=new double[cnum][fnum];
	    for(int i=0;i<cnum;i++) {
	    	for(int j=0;j<fnum;j++) {
	    		dist[i][j]=1/calDistanceSquare(cus.get(i),faC.get(j));
	    	}
	    }
	    
	    for(int i=0;i<cnum;i++) {
	    	IloLinearNumExpr expr4=ilcplex.linearNumExpr();
	    	for(int j=0;j<fnum;j++) {
	    		for(int k=0;k<K;k++) {
	    			expr4.addTerm(u[k]*dist[i][j], w[j][k]);
	    		}
	    	}
	    	ilcplex.addEq(b[i], ilcplex.prod(ilcplex.sum(t[i],expr4),1/t[i]));
	    }
	    
	    for(int i=0;i<cnum;i++) {
	    	ilcplex.addGe(ilcplex.prod(zeta[i], b[i]), 1); 
	    }
	      
	    ilcplex.addGe(omega, ilcplex.scalProd(h, zeta));
	    
	    //ilcplex.exportModel("SOCPmodel.lp");
	    
	    if(ilcplex.solve()) {
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
	    ilcplex.end();
	    return sol;
	}
	
	
//	public static void main(String[] args) throws IOException, IloException{
//		
//		Main m=new Main();
//		String pathname2="datasets/Facility_Leader.txt";
//		File filename2 = new File(pathname2);
//		String pathname3="datasets/Facility_Follower.txt";
//		File filename3 = new File(pathname3);
//		String pathname1="RandomData/"+"EX1/"+"30-20"+"/"+"30-20"+"-"+1+".txt";
//		File filename1 = new File(pathname1);
//		Instance instance=new Instance();
//		instance=m.initData(filename1,filename2,filename3);
//		
//		int fnum=instance.getFnum();
//		int cnum=instance.getCnum();
//		int R=2;
//		double V=200;
//		double[] h=new double[cnum];
//		for(int k=0;k<cnum;k++) {
//			h[k]=1.00/cnum;
//		}
//		double[] a= new double[fnum];
//		for(int j=0;j<fnum;j++) {
//			a[j]=50;
//		}
//		double[][] LeaderSol= {{0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0.0, 50.0, 0.0, 50.0, 0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0}};
//		double time1 = System.nanoTime();
//		
//		double beta =1;
//		double limit=0.1;
//		double[][] sol=SeparationSOCP(instance,LeaderSol,a,h,V,R,beta,limit);
//
//		//double[][] sol=Separation.SeparationProblem(LeaderSol, instance, V, a, h, R);
//		
//		double time2 = System.nanoTime();                                   
//		double time = (time2 - time1) / 1e9;//求解时间，单位s
//		
//		System.out.println(Arrays.toString(sol[0]));
//
//		System.out.println(Arrays.toString(sol[1]));
//		
//		System.out.println(time);
//		
//		
//		//System.out.println(Arrays.toString(bound));
//	}
}