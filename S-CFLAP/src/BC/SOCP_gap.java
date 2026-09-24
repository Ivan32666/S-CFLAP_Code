package BC;

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


public class SOCP_gap {
	
//	public static double calDistanceSquare(Customer cus,Facility fa) {
//		double ds=0;
//		ds=(cus.getX()-fa.getX())*(cus.getX()-fa.getX())+(cus.getY()-fa.getY())*(cus.getY()-fa.getY());	
//		return ds;
//	}
	public static double calDistanceSquare(Customer cus,Facility fa,double beta) {
		double ds=0;
		ds=Math.sqrt((cus.getX()-fa.getX())*(cus.getX()-fa.getX())+(cus.getY()-fa.getY())*(cus.getY()-fa.getY()));	
		double ds1=0;
		ds1=Math.pow(ds, beta);
		return ds1;
	}
	
	public static double[] calt(Instance instance, double[] u,double beta){
		int cnum=instance.getCnum();
		double[] t=new double[cnum];
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faC=instance.getFcandidate();
		for(int i=0;i<cnum;i++) {
			for(int k=0;k<faC.size();k++) {
				t[i]+=u[k]/calDistanceSquare(cus.get(i),faC.get(k),beta);
			}
		}
		
		return t;
	}
	
	

	public static double[][] SeparationSOCP(Instance instance,double[][] LeaderSol,double[] a,double[] h,double V,int R,double beta,double limit) throws IloException {
		IloCplex ilcplex =new IloCplex();
		ilcplex.setOut(null);
		ilcplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1.0e-1);
		int fnum=instance.getFnum();
		int cnum=instance.getCnum();
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faC=instance.getFcandidate();
		
	    IloIntVar[] y=new IloIntVar[fnum];  
	    IloNumVar[] v=new IloNumVar[fnum];
	    IloNumVar[] zeta=new IloNumVar[cnum];  
	    IloNumVar[] f1=new IloNumVar[cnum];
	    IloNumVar[] f2=new IloNumVar[cnum];  
	    IloNumVar[] b=new IloNumVar[cnum];
	    //IloIntVar[] tl=new IloIntVar[fnum];
	    
	    double[][] sol=new double[2][fnum];
	    for(int i=0;i<fnum;i++) {
	    	y[i]=ilcplex.boolVar("y"+i);
	    	v[i]=ilcplex.numVar(0, a[i],"v"+i);
	    	//tl[i]=ilcplex.boolVar();
	    }
	    for(int i=0;i<cnum;i++) {
	    	zeta[i]=ilcplex.numVar(0, 1,"zeta"+i);
	    	f1[i]=ilcplex.numVar(2, Double.MAX_VALUE,"f1"+i);
	    	f2[i]=ilcplex.numVar(0, Double.MAX_VALUE,"f2"+i);
	    	b[i]=ilcplex.numVar(1, Double.MAX_VALUE,"b"+i);
	    }

	  //objective
	    IloNumVar omega=ilcplex.numVar(0, 1);
	    
	    ilcplex.addMinimize(omega);
	  //constraints
	    
	    
	    IloLinearNumExpr expr2=ilcplex.linearNumExpr();
	    for(int i=0;i<fnum;i++) {
	    	expr2.addTerm(1, v[i]);
	    }
	    ilcplex.addLe(expr2, V);
	    
	    IloLinearNumExpr expr3=ilcplex.linearNumExpr();
	    for(int i=0;i<fnum;i++) {
	    	expr3.addTerm(1, y[i]);
	    }
	    ilcplex.addLe(expr3, R);
	   
	   
	    
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addLe(v[i],ilcplex.prod(a[i], y[i]));
	    }
	    
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addLe(y[i],1-LeaderSol[0][i]);
	    }
	    
	    double[] t=calt(instance,LeaderSol[1],beta);
	    
	    IloNumExpr[] expr4=new IloNumExpr[cnum];
	    
	    double[][] dist=new double[cnum][fnum];
	    for(int i=0;i<cnum;i++) {
	    	for(int j=0;j<fnum;j++) {
	    		dist[i][j]=1/calDistanceSquare(cus.get(i),faC.get(j),beta);
	    	}
	    }
	    
	    for(int i=0;i<cnum;i++) {
	    	expr4[i]=ilcplex.linearNumExpr();
	    	expr4[i]=ilcplex.scalProd(v, dist[i]);
	    	ilcplex.addEq(b[i], ilcplex.prod(ilcplex.sum(t[i],ilcplex.prod(1, expr4[i])),1/t[i]));
	    }
	    
	    for(int i=0;i<cnum;i++) {
	    	ilcplex.addEq(f1[i], ilcplex.sum(zeta[i],b[i]));
	    	ilcplex.addEq(f2[i], ilcplex.diff(b[i],zeta[i]));
	    }
	    
	    for(int i=0;i<cnum;i++) {
	    	ilcplex.addGe(ilcplex.prod(f1[i], f1[i]), ilcplex.sum(ilcplex.prod(f2[i], f2[i]),4));
	    }
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addGe(v[i], ilcplex.prod(limit*a[0], y[i]));
	    }
	    
	    ilcplex.addGe(omega, ilcplex.scalProd(h, zeta));
	    
	    //ilcplex.exportModel("SOCPmodel.lp");
	    
	    if(ilcplex.solve()) {
	    	double[] yVal = ilcplex.getValues(y);
	    	double[] vVal = ilcplex.getValues(v);
	    	sol[0]=yVal;
	    	//sol[1]=vVal;
	    	sol[1]=decimal(vVal);
	    	
	    }
	    ilcplex.end();
	    return sol;
	}
	
	public static double[] decimal(double[] vSol) {
		int num=vSol.length;
		double[] v=new double[num];
		for(int i=0;i<num;i++) {
			DecimalFormat df=new DecimalFormat("#.00000");
	        String format = df.format(vSol[i]);
	        v[i] = Double.parseDouble(format);
		}
		return v;
	}
	
	public static void main(String[] args) throws IOException, IloException{
		
		Main m=new Main();
		String pathname2="datasets/Facility_Leader.txt";
		File filename2 = new File(pathname2);
		String pathname3="datasets/Facility_Follower.txt";
		File filename3 = new File(pathname3);
		String pathname1="RandomData/"+"EX10/"+"100-20"+"/"+"100-20"+"-"+1+".txt";
		File filename1 = new File(pathname1);
		Instance instance=new Instance();
		instance=m.initData(filename1,filename2,filename3);
		
		int fnum=instance.getFnum();
		int cnum=instance.getCnum();
		int R=6;
		double V=200;
		double[] h=new double[cnum];
		for(int k=0;k<cnum;k++) {
			h[k]=1.00/cnum;
		}
		double[] a= new double[fnum];
		for(int j=0;j<fnum;j++) {
			a[j]=50;
		}
		double[][] LeaderSol= {{0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},{0.0, 20.0, 0.0, 20.0, 50.0, 0.0, 20.0, 0.0, 0.0, 0.0, 43.69103413317162, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 46.30896586682838}};
		double time1 = System.nanoTime();
		
		double beta =1;
		double limit=0.1;
		double[][] sol=SeparationSOCP(instance,LeaderSol,a,h,V,R,beta,limit);

		//double[][] sol=Separation.SeparationProblem(LeaderSol, instance, V, a, h, R);
		
		double time2 = System.nanoTime();                                   
		double time = (time2 - time1) / 1e9;//求解时间，单位s
		
		System.out.println(Arrays.toString(sol[0]));

		System.out.println(Arrays.toString(sol[1]));
		
		System.out.println(time);
		
		
		//System.out.println(Arrays.toString(bound));
	}
}