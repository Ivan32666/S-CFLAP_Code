package BC_Discrete;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import com.mathworks.toolbox.javabuilder.MWException;

import BC.SOCP;
import BC.UpAlg;
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

public class UpAlg_dis {
	
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
	
	public static double CalculateGradient(Instance instance,int j,int k,double[][] z,double[][] w,double[] u,double[] h) throws IloException {
		double gradient_j=0;
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faC=instance.getFcandidate();
		for(int i=0;i<cus.size();i++) {
			double sijk=-h[i]*u[k]/calDistanceSquare(cus.get(i),faC.get(j));
			double shareL=0;

			for(int m=0;m<faC.size();m++) {
				for(int n=0;n<u.length;n++) {
					shareL+=z[m][n]*u[n]/calDistanceSquare(cus.get(i),faC.get(m));
				}
				
			}
			double shareF=0;
			
			for(int m=0;m<faC.size();m++) {
				for(int n=0;n<u.length;n++) {
					shareF+=w[m][n]*u[n]/calDistanceSquare(cus.get(i),faC.get(m));
				}
				
			}
			
			gradient_j+=shareF*sijk/(shareL+shareF)/(shareL+shareF);

		}
		return gradient_j;
	}
//	public static double calDistanceSquare(Customer cus,Facility fa) {
//		double ds=0;
//		ds=Math.sqrt((cus.getX()-fa.getX())*(cus.getX()-fa.getX())+(cus.getY()-fa.getY())*(cus.getY()-fa.getY()));	
//		double ds1=0;
//		ds1=Math.pow(ds, 2);
//		return ds1;
//	}
	public static double calDistanceSquare(Customer cus,Facility fa) {
		double ds=0;
		ds=(cus.getX()-fa.getX())*(cus.getX()-fa.getX())+(cus.getY()-fa.getY())*(cus.getY()-fa.getY());	
		return ds;
	}
	
	
	public static double calMarketShare(double[] x,double[] y,Instance instance,double[][] z,double[][] w,double[] u,double[] h) {
		double share=0;
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faC=instance.getFcandidate();
		for(int i=0;i<cus.size();i++) {
			double shareL=0;

			for(int j=0;j<faC.size();j++) {
				for(int k=0;k<u.length;k++) {
					shareL+=z[j][k]*u[k]/calDistanceSquare(cus.get(i),faC.get(j));
				}
				
			}
			//System.out.println(shareL);
			double shareF=0;
			for(int j=0;j<faC.size();j++) {
				for(int k=0;k<u.length;k++) {
					shareF+=w[j][k]*u[k]/calDistanceSquare(cus.get(i),faC.get(j));
					
				}
			}
			
			//System.out.println(shareF);
			share+= h[i]*shareF/(shareL+shareF);
			
		}
		
		for(int j=0;j<faC.size();j++) {
			share-=x[j]*y[j];
		}

		return share;
	}
	public static Sol solveZeroLow(Instance instance,double[] x,double[][] z,double[] u,double[] a,int R,double V) throws IloException {
		int fnum=instance.getFnum();
		int cnum=instance.getCnum();
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faC=instance.getFcandidate();
		
		Sol sol =new Sol();
		IloCplex ilcplex =new IloCplex();
		ilcplex.setOut(null);
		
		int K=u.length;
	    IloIntVar[] y=new IloIntVar[fnum];  
	    IloIntVar[][] w=new IloIntVar[fnum][K];
	    //IloIntVar[] tl=new IloIntVar[fnum];  
	    
	    for(int i=0;i<fnum;i++) {
	    	y[i]=ilcplex.boolVar();
	    	for(int k=0;k<K;k++) {
	    		w[i][k]=ilcplex.boolVar();
	    	}
	    	
	    	//tl[i]=ilcplex.boolVar();
	    }
	    

	  //objective

	    IloNumExpr obj=ilcplex.numExpr();
	    
	    double[][] dist=new double[cnum][fnum];
	    for(int i=0;i<cnum;i++) {
	    	for(int j=0;j<fnum;j++) {
	    		dist[i][j]=1/calDistanceSquare(cus.get(i),faC.get(j));
	    	}
	    }
	    
	    for(int i=0;i<cnum;i++) {
	    	for(int j=0;j<fnum;j++) {
	    		for(int k=0;k<K;k++) {
	    			obj=ilcplex.sum(obj,ilcplex.prod(w[j][k],u[k]*dist[i][j]));
	    		}
	    	}
	    }
	  
	    ilcplex.addMaximize(obj);
	  //constraints
	    
	    IloLinearNumExpr expr2=ilcplex.linearNumExpr();
	    for(int i=0;i<fnum;i++) {
	    	expr2.addTerm(1, y[i]);
	    }
	    ilcplex.addLe(expr2, R);
	    
	    
	    IloLinearNumExpr expr3=ilcplex.linearNumExpr();
	    for(int j=0;j<fnum;j++) {
	    	for(int k=0;k<K;k++) {
	    		expr3.addTerm(w[j][k], u[k]);
	    	}
	    	
	    }
	    ilcplex.addLe(expr3, V);
	    

	    for(int i=0;i<fnum;i++) {
	    	ilcplex.addLe(y[i],1-x[i]);
	    }
	    
	    for(int j=0;j<fnum;j++) {
	    	IloLinearNumExpr expr4=ilcplex.linearNumExpr();
	    	for(int k=0;k<K;k++) {
		    	expr4.addTerm(1, w[j][k]);
		    }
	    	ilcplex.addEq(expr4,y[j]);
	    }

	    //ilcplex.exportModel("SP.lp");
	    if(ilcplex.solve()) {
	    	double[] yVal = ilcplex.getValues(y);
	    	double[][] wVal = new double[fnum][K];
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
	
	public static Cuts makecuts1(IloNumVar[] x, double[] xSol, IloModeler ilcplex,IloNumVar omega, Instance instance,double omega0,IloNumVar[][] z,double[][] zSol,double[] u,double[] h,double[] a,double V,int R) throws IloException {
		int fnum=instance.getFnum();
		Cuts cuts=new Cuts();
		ArrayList<IloNumExpr> cutLhs = new ArrayList<IloNumExpr>();
		ArrayList<Double> cutRhs = new ArrayList<Double>();
		
		IloLinearNumExpr lhsExpr=ilcplex.linearNumExpr();
		double rhsExpr=0;
		double lhs = 0, rhs = 0;
		ArrayList<Facility> faC=instance.getFcandidate();
		
		Sol upSol=new Sol();
		//double[] xSol=modifySol(xSol0,zSol);
		upSol.setLocation(xSol);
		upSol.setDesign_discrete(zSol);
		Sol lowSol=new Sol();
		
		if(isZero(zSol)==true) {
			lowSol=solveZeroLow(instance,xSol,zSol,u,a,R,V);
		}
		else {
			double time1=System.nanoTime();
			lowSol=SOCP_gap_dis.SeparationSOCP(instance, upSol, u, a, h, V, R);
			//lowSol=OA_gap.SeparationProblem(upSol, instance, V, u,a, h, R);
			double time2=System.nanoTime();
	 		time+=(time2 - time1) / 1e9;
		}
		lhsExpr.addTerm(1, omega); 
		
		rhsExpr+=calMarketShare(xSol,lowSol.getLocation(),instance,zSol,lowSol.getDesign_discrete(),u,h);
		
		for(int j=0;j<faC.size();j++) {
			for(int k=0;k<u.length;k++) {
				lhsExpr.addTerm(-CalculateGradient(instance,j,k,zSol,lowSol.getDesign_discrete(),u,h), z[j][k]);
				rhsExpr-=CalculateGradient(instance,j,k,zSol,lowSol.getDesign_discrete(),u,h)*zSol[j][k];
			}
			lhsExpr.addTerm(lowSol.getLocation()[j], x[j]);
			rhsExpr+=xSol[j]*lowSol.getLocation()[j];
		}
		lhs=omega0;
		rhs=calMarketShare(xSol,lowSol.getLocation(),instance,zSol,lowSol.getDesign_discrete(),u,h);
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
	
	public static Cuts makecuts2(IloNumVar[] x, double[] xSol, IloModeler ilcplex,IloNumVar omega, Instance instance,double omega0,IloNumVar[][] z,double[][] zSol,double[] u,double[] h,double[] a,double V,int R) throws IloException {
		int fnum=instance.getFnum();
		Cuts cuts=new Cuts();
		ArrayList<IloNumExpr> cutLhs = new ArrayList<IloNumExpr>();
		ArrayList<Double> cutRhs = new ArrayList<Double>();
		
		IloLinearNumExpr lhsExpr=ilcplex.linearNumExpr();
		double rhsExpr=0;
		double lhs = 0, rhs = 0;
		ArrayList<Facility> faC=instance.getFcandidate();
		
		Sol upSol=new Sol();
		//double[] xSol=modifySol(xSol0,zSol);
		upSol.setLocation(xSol);
		upSol.setDesign_discrete(zSol);
		Sol lowSol=new Sol();
		
		if(isZero(zSol)==true) {
			lowSol=solveZeroLow(instance,xSol,zSol,u,a,R,V);
		}
		else {
			double time1=System.nanoTime();
			lowSol=SOCP_dis.SeparationSOCP(instance, upSol,u, a, h, V, R);
			
			//lowSol=OA.SeparationProblem(upSol, instance, V, u,a, h, R);
			double time2=System.nanoTime();
	 		time+=(time2 - time1) / 1e9;
		}
		lhsExpr.addTerm(1, omega); 
		
		rhsExpr+=calMarketShare(xSol,lowSol.getLocation(),instance,zSol,lowSol.getDesign_discrete(),u,h);
		
		for(int j=0;j<faC.size();j++) {
			for(int k=0;k<u.length;k++) {
				lhsExpr.addTerm(-CalculateGradient(instance,j,k,zSol,lowSol.getDesign_discrete(),u,h), z[j][k]);
				rhsExpr-=CalculateGradient(instance,j,k,zSol,lowSol.getDesign_discrete(),u,h)*zSol[j][k];
			}
			lhsExpr.addTerm(lowSol.getLocation()[j], x[j]);
			rhsExpr+=xSol[j]*lowSol.getLocation()[j];
		}
		lhs=omega0;
		rhs=calMarketShare(xSol,lowSol.getLocation(),instance,zSol,lowSol.getDesign_discrete(),u,h);
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
	
	public static boolean isZero(double[][] z) {
		int num1=z.length;
		int num2=z[0].length;
		double sum=0;
		boolean flag;
		for(int i=0;i<num1;i++) {
			for(int j=0;j<num2;j++) {
				sum+=z[i][j];
			}
		}
		if(sum==0) {
			flag=true;
		}else {
			flag=false;
		}
		return flag;
	}

	public static Sol RMP(Instance instance,double U,double[] u,double[] a,double[] h,double V,int P,int R) throws IloException {
		IloCplex ilcplex =new IloCplex();
		ilcplex.setOut(null);
		//ilcplex.setParam(IloCplex.DoubleParam.TimeLimit, 10800);
		int fnum=instance.getFnum();
		int K=u.length;
	    IloIntVar[] x=new IloIntVar[fnum];  
	    IloIntVar[][] z=new IloIntVar[fnum][K];
	    //IloIntVar[] b=new IloIntVar[fnum];
	    
	    for(int i=0;i<fnum;i++) {
	    	x[i]=ilcplex.boolVar("x"+i);
	    	//b[i]=ilcplex.boolVar();
	    }
	    for(int j=0;j<fnum;j++) {
	    	for(int k=0;k<K;k++) {
	    		z[j][k]=ilcplex.boolVar("z"+j+k);
	    	}
	    }
	   
	  //objective
	    IloNumVar omega=ilcplex.numVar(0, 1);
	    
	    ilcplex.addMinimize(omega);
	  //constraints
	    
	    
	    IloLinearNumExpr expr2=ilcplex.linearNumExpr();
	    for(int j=0;j<fnum;j++) {
	    	for(int k=0;k<K;k++) {
	    		expr2.addTerm(u[k], z[j][k]);
	    	}
	    	
	    }
	    ilcplex.addLe(expr2, U);
	    
	    IloLinearNumExpr expr3=ilcplex.linearNumExpr();
	    for(int i=0;i<fnum;i++) {
	    	expr3.addTerm(1, x[i]);
	    }
	    ilcplex.addLe(expr3, P);
	    
	    for(int j=0;j<fnum;j++) {
	    	IloLinearNumExpr expr4=ilcplex.linearNumExpr();
	    	for(int k=0;k<K;k++) {
	    		expr4.addTerm(1, z[j][k]);
	    	}
	    	ilcplex.addEq(expr4, x[j]);
	    }
	    
	    
 		ilcplex.use(new LazyCallback(x,ilcplex,instance,omega,z,u,h,a,V,R));
 		
 		Sol sol=new Sol();
 		
	    if(ilcplex.solve()) {
	    	double[] xVal = ilcplex.getValues(x);
	    	double[][] zVal = new double[fnum][K];
	    	for(int j=0;j<fnum;j++) {
	    		for(int k=0;k<K;k++) {
	    			zVal[j][k]=ilcplex.getValue(z[j][k]);
	    		}
	    	}
	    	sol.setLocation(xVal);
	    	sol.setDesign_discrete(zVal);
	    
//	    	double gap=ilcplex.getMIPRelativeGap();
//	    	sol[2][0]=gap;
//	    	System.out.println(gap);
//	    	ilcplex.exportModel("Discretemodel.lp");
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
		IloIntVar[][] z;
		double[] h;
		double[] u;
		double[] a;
		double V;
		int R;
		
		LazyCallback(IloIntVar[] x0,IloCplex ilcplex0,Instance instance0,IloNumVar omega0,IloIntVar[][] z0,double[] u0,double[] h0,double[] a0,double V0,int R0){
			x=x0;
			ilcplex=ilcplex0;
			instance=instance0;
			omega=omega0;
			z=z0;
			u=u0;
			h=h0;
			a=a0;
			V=V0;
			R=R0;
		}
		
		public void main() throws IloException{
			cutNum++;
			double[] xSol =getValues(x);
			double omegaSol=getValue(omega);
			double[][] zSol=new double[xSol.length][u.length];
			for(int j=0;j<xSol.length;j++) {
				for(int k=0;k<u.length;k++) {
					zSol[j][k]=getValue(z[j][k]);
				}
			}
			
			
			cut1 = makecuts1(x, xSol,ilcplex,omega, instance,omegaSol,z,zSol,u,h,a,V,R);
			cutLhs = cut1.getLhs();
			cutRhs= cut1.getRhs();
			for(int i = 0; i< cutLhs.size(); i++) {
				addLocal(ilcplex.ge(ilcplex.diff(cutLhs.get(i),cutRhs.get(i)), 0));
			}
			if(cut1.getLhs().size()==0) {
				cut2 = makecuts2(x, xSol,ilcplex,omega, instance,omegaSol,z,zSol,u,h,a,V,R);
				cutLhs = cut2.getLhs();
				cutRhs= cut2.getRhs();
				for(int i = 0; i< cutLhs.size(); i++) {
					addLocal(ilcplex.ge(ilcplex.diff(cutLhs.get(i),cutRhs.get(i)), 0));
				}
			}
		}
	}
	public static double[] calDesign(double[] x,double[][] z, double[] u) {
		int fnum=x.length;
		int K=u.length;
		double[] design=new double[fnum];
		for(int j=0;j<fnum;j++) {
			if(Math.abs(x[j]-1)<1e-6) {
				for(int k=0;k<K;k++) {
					if(z[j][k]>0.1) {
						design[j]=u[k];	
					}
				}
			}
		}
		return design;
	}
	
	
//	public static double[] modifySol(double[] x,double[][] z) {
//		double[] newsol=x.clone();
//		
//		for(int j=0;j<x.length;j++) {
//			double sum=0;
//			for(int k=0;k<z[j].length;k++) {
//				sum+=z[j][k];
//			}
//			if(sum==0) {
//				newsol[j]=0;
//			}
//		}
//		return newsol;
//	}
	
	
	public static Solution solveOA(Instance instance,double[] h,double[] a,double[] u,int P,int R, double U, double V,double beta, double limit) throws IloException{
		Solution solution=new Solution();
		double time1 = System.nanoTime();
		Sol LeaderSol=RMP(instance,U,u,a,h,V,P,R);
//		for(int j=0;j<instance.getFnum();j++) {
//			System.out.println(Arrays.toString(LeaderSol.getDesign_discrete()[j]));
//		}
		double[][] LSol=new double[2][a.length];
		LSol[0]=LeaderSol.getLocation();
		LSol[1]=calDesign(LeaderSol.getLocation(),LeaderSol.getDesign_discrete(), u);
		
		//Sol FollwerSol=OA.SeparationProblem(LeaderSol, instance, V,u, a, h, R);
		SOCP sp=new SOCP();
		double[][] FollwerSol=sp.SeparationSOCP(instance, LSol, a, h, V, R,beta,limit);
		double time2 = System.nanoTime();                                   
		double Totaltime = (time2 - time1) / 1e9;//求解时间，单位s
		UpAlg mp=new UpAlg();
		double obj=mp.calMarketShare(LSol[0],FollwerSol[0],instance,LSol[1],FollwerSol[1],h,beta);
		solution.setLeaderLoca(LeaderSol.getLocation());
		solution.setLeaderAttrac(LSol[1]);
		solution.setFollowerLoca(FollwerSol[0]);
		solution.setFollowerAttrac(FollwerSol[1]);
		solution.setTime(Totaltime);
		solution.setObj(1-obj);
		solution.setCutnum(cutNum);
		cutNum=0;
		solution.setSPtime(time);
		time=0;
		return solution;
		
	}
	
	public static void main(String[] args) throws IOException, IloException, MWException{
		
		Main m=new Main();
		String pathname2="datasets/Facility_Leader.txt";
		File filename2 = new File(pathname2);
		String pathname3="datasets/Facility_Follower.txt";
		File filename3 = new File(pathname3);
		String pathname1="RandomData/EX1/20-10/20-10-1.txt";
		File filename1 = new File(pathname1);
		Instance instance=new Instance();
		instance=m.initData(filename1,filename2,filename3);
	
		int cnum=instance.getCnum();
		int fnum=instance.getCnum();
		int R=3;
		int P=3;
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
		double[] u= {20,30,40,50};
		Solution sol=UpAlg_dis.solveOA(instance,h,a,u,P,R, U, V,beta,limit);
		System.out.println(Arrays.toString(sol.getLeaderLoca()));
		System.out.println(Arrays.toString(sol.getLeaderAttrac()));
		System.out.println(Arrays.toString(sol.getFollowerLoca()));
		System.out.println(Arrays.toString(sol.getFollowerAttrac()));
		System.out.println(sol.getObj());
		System.out.println(sol.getTime());
		System.out.println(time);
		
		
		
//		double[][] LeaderSol= {{0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, -1.0, -0.0, 0.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, -0.0},{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 00.0, 50.0, 0.0, 0.0, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
//		double[][] FollowerSol= {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 50.0, 0, 49.99996, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
//		double beta=0.1;
//		double obj=calMarketShare(LeaderSol[0],FollowerSol[0],instance,LeaderSol[1],FollowerSol[1],h,beta);
	
		//double[][] sol=solveZeroLow(instance,LeaderSol[0],LeaderSol[1],a,R,V,beta);
//		System.out.println(1-obj);
		//System.out.println(Arrays.toString(sol[0]));
		//System.out.println(Arrays.toString(sol[1]));
		
//		[0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
//				[0.0, 20.0, 0.0, 20.0, 50.0, 0.0, 20.0, 0.0, 0.0, 0.0, 43.69103413317162, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 46.30896586682838]
//				[0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
//				[0.0, 0.0, 50.0, 0.0, 0.0, 50.0, 0.0, 0.0, 39.24206, 0.0, 0.0, 0.0, 20.75794, 20.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0]
//				0.5483196215823702
//				77.5968839
//				0.0
	}
 }
