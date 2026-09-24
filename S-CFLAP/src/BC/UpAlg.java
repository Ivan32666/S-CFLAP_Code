package BC;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import com.mathworks.toolbox.javabuilder.MWException;

import BC.OA.Callback;
import BC.OA.LazyCallback;
import Common.Customer;
import Common.Cuts;
import Common.Facility;
import Common.Instance;
import Common.Solution;
import Heuristic.Subproblem_VNS;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloModeler;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.LazyConstraintCallback;
import ilog.cplex.IloCplex.UserCutCallback;

public class UpAlg {
	
	public static int cutNum=0;
	public static double time=0;
	static Long seed;
    // 迭代次数
	static int epochs = 10;
    // 局部搜索次数
	static int localSearchCnt = 30;
    // VND参数
	static int lMax = 2;
    // VNS参数
	static int kMax = 2;
	
	public static double CalculateGradient(Instance instance,int j,double[] u,double[] v,double[] h,double beta) throws IloException {
		double gradient_j=0;
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faC=instance.getFcandidate();
		for(int i=0;i<cus.size();i++) {
			double wij=-h[i]/calDistanceSquare(cus.get(i),faC.get(j),beta);
			double shareL=0;

			for(int k=0;k<faC.size();k++) {
				shareL+=u[k]/calDistanceSquare(cus.get(i),faC.get(k),beta);
			}
			double shareF=0;

			for(int k=0;k<faC.size();k++) {
				shareF+=v[k]/calDistanceSquare(cus.get(i),faC.get(k),beta);
			}
			gradient_j+=shareF*wij/(shareL+shareF)/(shareL+shareF);

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
		ArrayList<Facility> faC=instance.getFcandidate();
		for(int i=0;i<cus.size();i++) {
			double shareL=0;

			for(int k=0;k<faC.size();k++) {
				shareL+=u[k]/calDistanceSquare(cus.get(i),faC.get(k),beta);
			}
			//System.out.println(shareL);
			double shareF=0;

			for(int k=0;k<faC.size();k++) {
				shareF+=v[k]/calDistanceSquare(cus.get(i),faC.get(k),beta);
			}
			//System.out.println(shareF);
			share+= h[i]*shareF/(shareL+shareF);
			
		}
		
		for(int j=0;j<faC.size();j++) {
			share-=x[j]*y[j];
		}

		return share;
	}
	public static double[][] solveZeroLow(Instance instance,double[] x,double[] u,double[] a,int R,double V,double beta,double limit) throws IloException {
		int fnum=instance.getFnum();
		int cnum=instance.getCnum();
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faC=instance.getFcandidate();
		
		double[][] sol =new double[2][fnum];
		IloCplex ilcplex =new IloCplex();
		ilcplex.setOut(null);
		
	    IloIntVar[] y=new IloIntVar[fnum];  
	    IloNumVar[] v=new IloNumVar[fnum];
	    //IloIntVar[] tl=new IloIntVar[fnum];  
	    
	    for(int i=0;i<fnum;i++) {
	    	y[i]=ilcplex.boolVar();
	    	v[i]=ilcplex.numVar(0, a[i]);
	    	//tl[i]=ilcplex.boolVar();
	    }
	    

	  //objective
	    IloNumExpr[] expr1=new IloNumExpr[cnum];
	    IloNumExpr obj=ilcplex.numExpr();
	    
	    double[][] dist=new double[cnum][fnum];
	    for(int i=0;i<cnum;i++) {
	    	for(int j=0;j<fnum;j++) {
	    		dist[i][j]=1/calDistanceSquare(cus.get(i),faC.get(j),beta);
	    	}
	    }
	    for(int i=0;i<cnum;i++) {
	    	expr1[i]=ilcplex.scalProd(v, dist[i]);
	    }
	    for(int j=0;j<cnum;j++) {
	    	obj=ilcplex.sum(obj,ilcplex.prod(1, expr1[j]));
	    }
	    
	    ilcplex.addMaximize(obj);
	  //constraints
	    
	    IloLinearNumExpr expr2=ilcplex.linearNumExpr();
	    for(int i=0;i<fnum;i++) {
	    	expr2.addTerm(1, y[i]);
	    }
	    ilcplex.addLe(expr2, R);
	    
	    
	    IloLinearNumExpr expr3=ilcplex.linearNumExpr();
	    for(int i=0;i<fnum;i++) {
	    	expr3.addTerm(1, v[i]);
	    }
	    ilcplex.addLe(expr3, V);
	    
 
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addLe(v[i],ilcplex.prod(a[i], y[i]));
	    }

	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addLe(y[i],1-x[i]);
	    }
	    
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addGe(v[i],ilcplex.prod(limit*a[0], y[i]));
	    }
	   
	    
	    //ilcplex.exportModel("SP.lp");
	    if(ilcplex.solve()) {
	    	double[] yVal = ilcplex.getValues(y);
	    	double[] vVal = ilcplex.getValues(v);
	    	sol[0]=yVal;
	    	sol[1]=vVal;
	    }
	    return sol;
	}
	
	public static Cuts makecuts1(IloNumVar[] x, double[] xSol, IloModeler ilcplex,IloNumVar omega, Instance instance,double omega0,IloNumVar[] u,double[] uSol,double[] h,double[] a,double V,int R,double beta,double limit) throws IloException {
		int fnum=instance.getFnum();
		Cuts cuts=new Cuts();
		ArrayList<IloNumExpr> cutLhs = new ArrayList<IloNumExpr>();
		ArrayList<Double> cutRhs = new ArrayList<Double>();
		
		IloLinearNumExpr lhsExpr=ilcplex.linearNumExpr();
		double rhsExpr=0;
		double lhs = 0, rhs = 0;
		ArrayList<Facility> faC=instance.getFcandidate();
		double[][] upSol=new double[2][fnum];
		upSol[0]=xSol;
		upSol[1]=uSol;
		double[][] lowSol=new double[2][fnum];
		if(isZero(uSol)==true) {
			lowSol=solveZeroLow(instance,xSol,uSol,a,R,V,beta,limit);
			 
		}
		else {
			double time1=System.nanoTime();
			lowSol=SOCP_gap.SeparationSOCP(instance, upSol, a, h, V, R,beta,limit);
			//lowSol=OA_gap.SeparationProblem(upSol, instance, V, a, h, R, limit);
//			Subproblem_VNS vns=new Subproblem_VNS(seed, lMax, kMax, epochs, localSearchCnt);
//			lowSol=vns.Separation_VNS(instance,R,xSol,uSol,h,V,a,limit);
			double time2=System.nanoTime();
	 		time+=(time2 - time1) / 1e9;
		}
		lhsExpr.addTerm(1, omega); 
		
		rhsExpr+=calMarketShare(xSol,lowSol[0],instance,uSol,lowSol[1],h,beta);
		
		
		for(int i=0;i<faC.size();i++) {
			lhsExpr.addTerm(-CalculateGradient(instance,i,uSol,lowSol[1],h,beta), u[i]);
			lhsExpr.addTerm(lowSol[0][i], x[i]);
			rhsExpr-=CalculateGradient(instance,i,uSol,lowSol[1],h,beta)*uSol[i];
			rhsExpr+=xSol[i]*lowSol[0][i];
		}
		lhs=omega0;
		rhs=calMarketShare(xSol,lowSol[0],instance,uSol,lowSol[1],h,beta);
//		System.out.println(lhs);
//		System.out.println(rhs);
		if(lhs<rhs&&(rhs-lhs>1e-6)) {
			cutLhs.add(lhsExpr);
			cutRhs.add(rhsExpr);
		}

		cuts.setLhs(cutLhs);
		cuts.setRhs(cutRhs);
		
		return cuts;
	}
	
	public static Cuts makecuts2(IloNumVar[] x, double[] xSol, IloModeler ilcplex,IloNumVar omega, Instance instance,double omega0,IloNumVar[] u,double[] uSol,double[] h,double[] a,double V,int R,double beta,double limit) throws IloException {
		int fnum=instance.getFnum();
		Cuts cuts=new Cuts();
		ArrayList<IloNumExpr> cutLhs = new ArrayList<IloNumExpr>();
		ArrayList<Double> cutRhs = new ArrayList<Double>();
		
		IloLinearNumExpr lhsExpr=ilcplex.linearNumExpr();
		double rhsExpr=0;
		double lhs = 0, rhs = 0;
		ArrayList<Facility> faC=instance.getFcandidate();
		double[][] upSol=new double[2][fnum];
		upSol[0]=xSol;
		upSol[1]=uSol;
		double[][] lowSol=new double[2][fnum];
		if(isZero(uSol)==true) {
			lowSol=solveZeroLow(instance,xSol,uSol,a,R,V,beta,limit);
			 
		}
		else {
			double time1=System.nanoTime();
			lowSol=SOCP.SeparationSOCP(instance, upSol, a, h, V, R,beta,limit);
		    //lowSol=OA.SeparationProblem(upSol, instance, V, a, h, R, limit);
//			Subproblem_VNS vns=new Subproblem_VNS(seed, lMax, kMax, epochs, localSearchCnt);
//			lowSol=vns.Separation_VNS(instance,R,xSol,uSol,h,V,a,limit);
			double time2=System.nanoTime();
	 		time+=(time2 - time1) / 1e9;
		}
		
		
		lhsExpr.addTerm(1, omega); 
		
		rhsExpr+=calMarketShare(xSol,lowSol[0],instance,uSol,lowSol[1],h,beta);
		
		
		for(int i=0;i<faC.size();i++) {
			lhsExpr.addTerm(-CalculateGradient(instance,i,uSol,lowSol[1],h,beta), u[i]);
			lhsExpr.addTerm(lowSol[0][i], x[i]);
			rhsExpr-=CalculateGradient(instance,i,uSol,lowSol[1],h,beta)*uSol[i];
			rhsExpr+=xSol[i]*lowSol[0][i];
		}
		lhs=omega0;
		rhs=calMarketShare(xSol,lowSol[0],instance,uSol,lowSol[1],h,beta);
//		System.out.println(lhs);
//		System.out.println(rhs);
		if(lhs<rhs&&(rhs-lhs>1e-6)) {
			cutLhs.add(lhsExpr);
			cutRhs.add(rhsExpr);
		}

		cuts.setLhs(cutLhs);
		cuts.setRhs(cutRhs);
		
		return cuts;
	}
	
	public static boolean isZero(double[] u) {
		int num=u.length;
		double sum=0;
		boolean flag;
		for(int i=0;i<num;i++) {
			sum+=u[i];
		}
		if(sum==0) {
			flag=true;
		}else {
			flag=false;
		}
		return flag;
	}

	public static double[][] RMP(Instance instance,double U,double[] a,double[] h,double V,int P,int R, double beta,double limit) throws IloException {
		IloCplex ilcplex =new IloCplex();
		ilcplex.setOut(null);
		//ilcplex.setParam(IloCplex.DoubleParam.TimeLimit, 10800);
		int fnum=instance.getFnum();
	    IloIntVar[] x=new IloIntVar[fnum];  
	    IloNumVar[] u=new IloNumVar[fnum];
	    //IloIntVar[] b=new IloIntVar[fnum];
	    
	    double[][] sol=new double[3][fnum];
	    for(int i=0;i<fnum;i++) {
	    	x[i]=ilcplex.boolVar();
	    	u[i]=ilcplex.numVar(0, a[i]);
	    	//b[i]=ilcplex.boolVar();
	    }
	   
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
	    ilcplex.addLe(expr3, P);
	    
	    
	    
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addLe(u[i],ilcplex.prod(a[i], x[i]));
	    }
	    
	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addGe(u[i],ilcplex.prod(limit*a[0], x[i]));
	    }
	    
 		ilcplex.use(new LazyCallback(x,ilcplex,instance,omega,u,h,a,V,R,beta,limit));
 		
	    if(ilcplex.solve()) {
	    	double[] xVal = ilcplex.getValues(x);
	    	double[] uVal = ilcplex.getValues(u);
	    	sol[0]=xVal;
	    	sol[1]=uVal;
	    	double gap=ilcplex.getMIPRelativeGap();
	    	sol[2][0]=gap;
	    	System.out.println(gap);
	    	//ilcplex.exportModel("model.lp");
	    }
	    ilcplex.end();
	    return sol;
	}
	
	public static class LazyCallback extends LazyConstraintCallback{
		Cuts cut1;
		Cuts cut2;
		ArrayList<IloNumExpr> cutLhs;
		ArrayList<Double> cutRhs;
		IloIntVar[] x;
		IloCplex ilcplex;
		Instance instance;
		IloNumVar omega;
		IloNumVar[] u;
		double[] h;
		double[] a;
		double V;
		int R;
		double beta;
		double limit;
		
		LazyCallback(IloIntVar[] x0,IloCplex ilcplex0,Instance instance0,IloNumVar omega0,IloNumVar[] u0,double[] h0,double[] a0,double V0,int R0,double beta0,double limit0){
			x=x0;
			ilcplex=ilcplex0;
			instance=instance0;
			omega=omega0;
			u=u0;
			h=h0;
			a=a0;
			V=V0;
			R=R0;
			beta=beta0;
			limit=limit0;
		}
		
		public void main() throws IloException{
			cutNum++;
			double[] xSol =getValues(x);
			double omegaSol=getValue(omega);
			double[] uSol=getValues(u);
			
			cut1 = makecuts1(x, xSol,ilcplex,omega, instance,omegaSol,u,uSol,h,a,V,R,beta,limit);
			
			cutLhs = cut1.getLhs();
			cutRhs= cut1.getRhs();
			for(int i = 0; i< cutLhs.size(); i++) {
				addLocal(ilcplex.ge(ilcplex.diff(cutLhs.get(i),cutRhs.get(i)), 0));
			}
			if(cut1.getLhs().size()==0) {
				cut2 = makecuts2(x, xSol,ilcplex,omega, instance,omegaSol,u,uSol,h,a,V,R,beta,limit);
				cutLhs = cut2.getLhs();
				cutRhs= cut2.getRhs();
				for(int i = 0; i< cutLhs.size(); i++) {
					addLocal(ilcplex.ge(ilcplex.diff(cutLhs.get(i),cutRhs.get(i)), 0));
				}
			}
		}
	}
	
	
	public static Solution solveOA(Instance instance,double[] h,double[] a,int P,int R, double U, double V,double beta,double limit) throws IloException{
		Solution solution=new Solution();
		double time1 = System.nanoTime();
		double[][] LeaderSol=RMP(instance,U,a,h,V,P,R,beta,limit);
		//double[][] FollwerSol=OA.SeparationProblem(LeaderSol, instance, V, a, h, R,limit);
		double[][] FollwerSol=SOCP.SeparationSOCP(instance, LeaderSol, a, h, V, R, beta,limit);
		double time2 = System.nanoTime();                                   
		double Totaltime = (time2 - time1) / 1e9;//求解时间，单位s
		double obj=calMarketShare(LeaderSol[0],FollwerSol[0],instance,LeaderSol[1],FollwerSol[1],h,beta);
		solution.setLeaderLoca(LeaderSol[0]);
		solution.setLeaderAttrac(LeaderSol[1]);
		solution.setFollowerLoca(FollwerSol[0]);
		solution.setFollowerAttrac(FollwerSol[1]);
		solution.setTime(Totaltime);
		solution.setObj(1-obj);
		solution.setCutnum(cutNum);
		cutNum=0;
		solution.setSPtime(time);
		time=0;
		solution.setGap(LeaderSol[2][0]);
		return solution;
		
	}
	
	public static void main(String[] args) throws IOException, IloException, MWException{
		
		Main m=new Main();
		String pathname2="datasets/Facility_Leader.txt";
		File filename2 = new File(pathname2);
		String pathname3="datasets/Facility_Follower.txt";
		File filename3 = new File(pathname3);
		String pathname1="RandomData/EX for Table 8/100-20-1.txt";
		File filename1 = new File(pathname1);
		Instance instance=new Instance();
		instance=m.initData(filename1,filename2,filename3);
	
		int cnum=instance.getCnum();
		int fnum=instance.getCnum();
		int R=2;
		int P=2;
		double V=200;
		double U=200;
		double[] h=new double[cnum];
		for(int k=0;k<cnum;k++) {
			h[k]=1.00/cnum;
		}
		double[] a= new double[fnum];
		for(int j=0;j<fnum;j++) {
			a[j]=100;
		}
		double limit=0.4;
		double beta=2;
		Solution sol=UpAlg.solveOA(instance,h,a,P,R,U,V,beta,limit);
		System.out.println(Arrays.toString(sol.getLeaderLoca()));
		System.out.println(Arrays.toString(sol.getLeaderAttrac()));
		System.out.println(Arrays.toString(sol.getFollowerLoca()));
		System.out.println(Arrays.toString(sol.getFollowerAttrac()));
		System.out.println(sol.getObj());
		System.out.println(sol.getTime());
		System.out.println(time);
		
		
		
	}
 }
