package BC;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;

import com.mathworks.toolbox.javabuilder.MWException;

import Common.Customer;
import Common.Facility;
import Common.Instance;
import Common.Solution;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

public class simultaneousCFLDP {
	
	public static double calDistanceSquare(Customer cus,Facility fa) {
		double ds=0;
		ds=(cus.getX()-fa.getX())*(cus.getX()-fa.getX())+(cus.getY()-fa.getY())*(cus.getY()-fa.getY());	
		return ds;
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
	
	public static double[][] calUtility(Instance instance, double[] u,double[] v){
		int fnum=instance.getFnum();
		int cnum=instance.getCnum();
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faC=instance.getFcandidate();
		double[][] utility=new double[2][cnum];
		for(int i=0;i<cnum;i++) {
			for(int j=0;j<fnum;j++) {
				utility[0][i]=utility[0][i]+u[j]/calDistanceSquare(cus.get(i),faC.get(j));
				utility[1][i]=utility[1][i]+v[j]/calDistanceSquare(cus.get(i),faC.get(j));
			}
		}
		return utility;
	}
	
	public static double[][] Up0(Instance instance,double[] a,double[] h,double U,double limit) throws IloException{
		IloCplex ilcplex =new IloCplex();
		ilcplex.setOut(null);
		int fnum=instance.getFnum();
		int cnum=instance.getCnum();
		double[][] leaderSol=new double[2][fnum];
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faC=instance.getFcandidate();
		
	    IloIntVar[] x=new IloIntVar[fnum];  
	    IloNumVar[] u=new IloNumVar[fnum];

	    for(int i=0;i<fnum;i++) {
	    	x[i]=ilcplex.boolVar("x"+i);
	    	u[i]=ilcplex.numVar(0, a[i],"u"+i);
	    	//tl[i]=ilcplex.boolVar();
	    }
	   
	    
	  //objective
	    IloLinearNumExpr obj=ilcplex.linearNumExpr();
	    for(int i=0;i<cnum;i++) {
	    	for(int j=0;j<fnum;j++) {
	    		obj.addTerm(h[i]/calDistanceSquare(cus.get(i),faC.get(j)), u[j]);
	    	}
	    }
	    
	    ilcplex.addMaximize(obj);
	  //constraints
	    
	    
	    IloLinearNumExpr expr2=ilcplex.linearNumExpr();
	    for(int i=0;i<fnum;i++) {
	    	expr2.addTerm(1, u[i]);
	    }
	    ilcplex.addLe(expr2, U);
	    
	    IloLinearNumExpr expr3=ilcplex.linearNumExpr();
	    for(int i=0;i<fnum;i++) {
	    	expr3.addTerm(1, x[i]);
	    }
	    ilcplex.addLe(expr3, 1);
	   
	   
	    
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addLe(u[i],ilcplex.prod(a[i], x[i]));
	    }
	    
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addGe(u[i], ilcplex.prod(limit*a[0], x[i]));
	    }
	    
	   
	    //ilcplex.exportModel("SOCPmodel.lp");
	    
	    if(ilcplex.solve()) {
	    	double[] xVal = ilcplex.getValues(x);
	    	double[] uVal = ilcplex.getValues(u);
	    	leaderSol[0]=xVal;
	    	leaderSol[1]=decimal(uVal);
	    	
	    }
	    ilcplex.end();
	    return leaderSol;
	}
	
	public static double[][] UP(Instance instance,double[][] Lsol,double[][] Fsol, double[] a,double[] h,double U,double limit) throws IloException{
		IloCplex ilcplex =new IloCplex();
		ilcplex.setOut(null);
		int fnum=instance.getFnum();
		int cnum=instance.getCnum();
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faC=instance.getFcandidate();
		
	    IloIntVar[] x=new IloIntVar[fnum];  
	    IloNumVar[] u=new IloNumVar[fnum];
	    IloNumVar[] zeta=new IloNumVar[cnum];  
	    IloNumVar[] f1=new IloNumVar[cnum];
	    IloNumVar[] f2=new IloNumVar[cnum];  
	    IloNumVar[] b=new IloNumVar[cnum];
	    
	    double[][] leaderSol=new double[2][fnum];
	    
	    for(int i=0;i<fnum;i++) {
	    	x[i]=ilcplex.boolVar("x"+i);
	    	u[i]=ilcplex.numVar(0, a[i],"u"+i);
	    	//tl[i]=ilcplex.boolVar();
	    }
	    for(int i=0;i<cnum;i++) {
	    	zeta[i]=ilcplex.numVar(0, 1,"zeta"+i);
	    	f1[i]=ilcplex.numVar(2, Double.MAX_VALUE,"f1"+i);
	    	f2[i]=ilcplex.numVar(0, Double.MAX_VALUE,"f2"+i);
	    	b[i]=ilcplex.numVar(1, Double.MAX_VALUE,"b"+i);
	    }
	    
	    double[][] utility=calUtility(instance, Lsol[1],Fsol[1]);
	    
	  //objective
	    IloNumVar omega=ilcplex.numVar(0, 1);
	    
	    ilcplex.addMinimize(omega);
	  //constraints
	    
	    
	    IloLinearNumExpr expr2=ilcplex.linearNumExpr();
	    for(int i=0;i<fnum;i++) {
	    	expr2.addTerm(1, u[i]);
	    }
	    ilcplex.addLe(expr2, U);
	    
	    IloLinearNumExpr expr3=ilcplex.linearNumExpr();
	    for(int i=0;i<fnum;i++) {
	    	expr3.addTerm(1, x[i]);
	    }
	    ilcplex.addLe(expr3, 1);
	   
	   
	    
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addLe(u[i],ilcplex.prod(a[i], x[i]));
	    }
	    
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addLe(x[i],1-Lsol[0][i]);
	    	ilcplex.addLe(x[i],1-Fsol[0][i]);
	    }
	    
	    double[] t=utility[1];
	    
	    IloNumExpr[] expr4=new IloNumExpr[cnum];
	    
	    double[][] dist=new double[cnum][fnum];
	    for(int i=0;i<cnum;i++) {
	    	for(int j=0;j<fnum;j++) {
	    		dist[i][j]=1/calDistanceSquare(cus.get(i),faC.get(j));
	    	}
	    }
	    
	    for(int i=0;i<cnum;i++) {
	    	expr4[i]=ilcplex.linearNumExpr();
	    	expr4[i]=ilcplex.scalProd(u, dist[i]);
	    	ilcplex.addEq(b[i], ilcplex.prod(ilcplex.sum(ilcplex.sum(t[i],expr4[i]),utility[1][i]),1/t[i]));
	    }
	    
	    for(int i=0;i<cnum;i++) {
	    	ilcplex.addEq(f1[i], ilcplex.sum(zeta[i],b[i]));
	    	ilcplex.addEq(f2[i], ilcplex.diff(b[i],zeta[i]));
	    }
	    
	    for(int i=0;i<cnum;i++) {
	    	ilcplex.addGe(ilcplex.prod(f1[i], f1[i]), ilcplex.sum(ilcplex.prod(f2[i], f2[i]),4));
	    }
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addGe(u[i], ilcplex.prod(limit*a[0], x[i]));
	    }
	    
	    ilcplex.addGe(omega, ilcplex.scalProd(h, zeta));
	    
	    //ilcplex.exportModel("SOCPmodel.lp");
	    
	    if(ilcplex.solve()) {
	    	double[] xVal = ilcplex.getValues(x);
	    	double[] uVal = ilcplex.getValues(u);
	    	leaderSol[0]=xVal;
	    	leaderSol[1]=decimal(uVal);
	    	
	    }
	    ilcplex.end();
	    return leaderSol;
	}
	
	public static double[][] SP(Instance instance,double[][] Lsol,double[][] Fsol, double[] a,double[] h,double V,double limit) throws IloException{
		IloCplex ilcplex =new IloCplex();
		ilcplex.setOut(null);
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
	    
	    double[][] followerSol=new double[2][fnum];
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
	    
	    double[][] utility=calUtility(instance, Lsol[1],Fsol[1]);
	    
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
	    ilcplex.addLe(expr3, 1);
	   
	   
	    
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addLe(v[i],ilcplex.prod(a[i], y[i]));
	    }
	    
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addLe(y[i],1-Lsol[0][i]);
	    	ilcplex.addLe(y[i],1-Fsol[0][i]);
	    }
	    
	    double[] t=utility[0];
	    
	    IloNumExpr[] expr4=new IloNumExpr[cnum];
	    
	    double[][] dist=new double[cnum][fnum];
	    for(int i=0;i<cnum;i++) {
	    	for(int j=0;j<fnum;j++) {
	    		dist[i][j]=1/calDistanceSquare(cus.get(i),faC.get(j));
	    	}
	    }
	    
	    for(int i=0;i<cnum;i++) {
	    	expr4[i]=ilcplex.linearNumExpr();
	    	expr4[i]=ilcplex.scalProd(v, dist[i]);
	    	ilcplex.addEq(b[i], ilcplex.prod(ilcplex.sum(ilcplex.sum(t[i],expr4[i]),utility[1][i]),1/t[i]));
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
	    	followerSol[0]=yVal;
	    	followerSol[1]=decimal(vVal);
	    	
	    }
	    ilcplex.end();
	    return followerSol;
	}
	
	public static double[][] combine(double[][] sol1, double[][] sol2){
		double[][] newsol=new double[2][sol1[1].length];
		for(int j=0;j<sol1[1].length;j++) {
			newsol[0][j]=sol1[0][j]+sol2[0][j];
			newsol[1][j]=sol1[1][j]+sol2[1][j];
		}
		return newsol;
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
			//System.out.println(shareL);
			double shareF=0;

			for(int k=0;k<faC.size();k++) {
				shareF+=v[k]/calDistanceSquare(cus.get(i),faC.get(k));
			}
			//System.out.println(shareF);
			share+= h[i]*shareF/(shareL+shareF);
			
		}

		return share;
	}
	
	public static Solution simultaneousModel(Instance instance,double[] h,double[] a,int P,int R, double U, double V,double beta,double limit) throws IloException {
		Solution solution=new Solution();
		double time1 = System.nanoTime();
		int fnum=instance.getFnum();
		double[][] LeaderSol1=Up0(instance,a,h,U/2,limit);
		double[][] Fsol= new double[2][fnum];
		double[][] FollwerSol1=SP(instance,LeaderSol1,Fsol, a,h,V/2,limit);
		double[][] LeaderSol2=UP(instance,LeaderSol1,FollwerSol1, a, h, U/2,limit);
		double[][] FollwerSol2=SP(instance,combine(LeaderSol1,LeaderSol2),FollwerSol1, a,h,V/2,limit);
		
		double[][] LeaderSol=combine(LeaderSol1, LeaderSol2);
		double[][] FollowerSol=combine(FollwerSol1, FollwerSol2);
		//double[][] FollowerSol=SOCP.SeparationSOCP(instance, LeaderSol, a, h, V, R, beta,limit);
		double time2 = System.nanoTime();                                   
		double Totaltime = (time2 - time1) / 1e9;//求解时间，单位s
		double obj=calMarketShare(LeaderSol[0],FollowerSol[0],instance,LeaderSol[1],FollowerSol[1],h);
		solution.setLeaderLoca(LeaderSol[0]);
		solution.setLeaderAttrac(LeaderSol[1]);
		solution.setFollowerLoca(FollowerSol[0]);
		solution.setFollowerAttrac(FollowerSol[1]);
		solution.setTime(Totaltime);
		solution.setObj(1-obj);
		return solution;
	}
	
	public static void main(String[] args) throws IOException, IloException, MWException{
		
		Main m=new Main();
		String pathname2="datasets/Facility_Leader.txt";
		File filename2 = new File(pathname2);
		String pathname3="datasets/Facility_Follower.txt";
		File filename3 = new File(pathname3);
		String pathname1="RandomData/EX for Table 8/200-20-1.txt";
		File filename1 = new File(pathname1);
		Instance instance=new Instance();
		instance=m.initData(filename1,filename2,filename3);
	
		int cnum=instance.getCnum();
		int fnum=instance.getCnum();
		int R=2;
		int P=2;
		double V=100;
		double U=100;
		double[] h=new double[cnum];
		for(int k=0;k<cnum;k++) {
			h[k]=1.00/cnum;
		}
		double[] a= new double[fnum];
		for(int j=0;j<fnum;j++) {
			a[j]=50;
		}
		double limit=0.4;
		double beta=1;
		Solution sol=simultaneousModel(instance,h,a,P,R, U, V,beta,limit);
		System.out.println(Arrays.toString(sol.getLeaderLoca()));
		System.out.println(Arrays.toString(sol.getLeaderAttrac()));
		System.out.println(Arrays.toString(sol.getFollowerLoca()));
		System.out.println(Arrays.toString(sol.getFollowerAttrac()));
		System.out.println(sol.getObj());
		System.out.println(sol.getTime());
		

	}
}
