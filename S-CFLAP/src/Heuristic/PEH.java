package Heuristic;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import BC.Main;
import BC.SOCP;
import Common.Customer;
import Common.Facility;
import Common.Instance;
import Common.Solution;
import ilog.concert.IloException;

public class PEH {
	 public static int interation=0;
    
	 private List<List<Integer>> resc = new ArrayList<>();
		
	 public List<List<Integer>> combine(int n, int k) {
		
		 if (n <= 0 || k <= 0 || k > n) {
			return resc;
		}
		
		List<Integer> c = new ArrayList<>();
		generateCombinations(n, k, 1, c);
		return resc;
	
		}

	    /**
		 * 回溯求所有组合结果
		 * @param n 
		 * @param k 
		 * @param start 开始搜索新元素的位置
		 * @param c 当前已经找到的组合
		 */
		private void generateCombinations(int n,int k,int start,List<Integer> c) {
			if (c.size() == k) {
			    //这里需要注意java的值传递
				//此处必须使用重新创建对象的形式，否则 res 列表中存放的都是同一个引用
				resc.add(new ArrayList<>(c));
				return;
			}
			
			//通过终止条件，进行剪枝优化，避免无效的递归
			//c中还剩 k - c.size()个空位，所以[ i ... n]中至少要有k-c.size()个元素
			//所以i最多为 n - (k - c.size()) + 1
			for(int i = start;i <= n - (k - c.size()) + 1; i++) {
				c.add(i);
				generateCombinations(n, k, i + 1, c);
				//记得回溯状态啊
				c.remove(c.size() - 1);
			}
		}
	
		public static int[] enumrate(int fnum,int p,List<Integer> l) {
			int[] x=new int[fnum];
			for(int i=0;i<fnum;i++) {
				for(int j=0;j<l.size();j++) {
					if(l.get(j)==i+1) {
						x[i]=1;
					}
				}
			}
			return x;
		}
		public static FollowerSol subproblem(Instance instance,int[] x,double[] u,int P,double U,double[] a,double[] h, double limit) {
			PEH p=new PEH();
			int fnum=instance.getFnum();
			int[] loc=new int[fnum];
			double[] attrac=new double[fnum];
			double marketshare=1;
			FollowerSol fsol=new FollowerSol(marketshare,loc,attrac);
			List<List<Integer>> l=p.combine(fnum, P);
			 for(int i=0;i<l.size();i++) {
				 interation++;
//				 System.out.println(interation);
				 int[] y=enumrate(fnum,P,l.get(i));
				 if(isValid(x,y)==false) {
					 continue;
				 }
				 double[] v=null;
			     PSO pso=new PSO();
				 v=pso.U_sol(instance, y, U, a,limit);
				 //System.out.println(Arrays.toString(v));
				 marketshare=calMarketShare(u,v,instance,h);
				 if(marketshare<fsol.getMarketshare()) {
					 fsol.setLocation(y);
					 fsol.setAttract(v);
					 fsol.setMarketshare(marketshare);
				 }
			 }
			return fsol;
		}
		public static double calDistanceSquare(Customer cus,Facility fa) {
			double ds=0;
			ds=(cus.getX()-fa.getX())*(cus.getX()-fa.getX())+(cus.getY()-fa.getY())*(cus.getY()-fa.getY());	
			return ds;
		}
		public static double calMarketShare(double[] u,double[] v,Instance instance,double[] h) {
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
				for(int k=0;k<faC.size();k++) {
					shareL+=u[k]/calDistanceSquare(cus.get(i),faC.get(k));
				}
				double shareF=0;
				for(int k=0;k<faF0.size();k++) {
					shareF+=faF0.get(k).getSize()/calDistanceSquare(cus.get(i),faF0.get(k));
				}
				for(int k=0;k<faC.size();k++) {
					shareF+=v[k]/calDistanceSquare(cus.get(i),faC.get(k));
				}
				share+=h[i]*shareL/(shareL+shareF);
			}
			return share;
		}
		public static boolean isValid(int[] s1,int[] s2) {
	    	boolean flag=true;
	    	for(int i=0;i<s1.length;i++) {
	    		if(s1[i]==1&&s2[i]==1) {
	    			flag=false;
	    			break;
	    		}
	    	}
	    	return flag;
	    }
		public static Sol PEHAlg(Instance instance,int P,int R,double U,double V,double[] a,double[] h, double limit) {
			PEH p=new PEH();
			int fnum=instance.getFnum();
			int[] leaderloc=new int[fnum];
			int[] followerloc=new int[fnum];
			double[] leaderattrac=new double[fnum];
			double[] followerattrac=new double[fnum];
			double marketshare=0;
			Sol lsol=new Sol(marketshare,leaderloc,followerloc,leaderattrac,followerattrac);
			List<List<Integer>> l=p.combine(fnum, P);
			double time1 = System.nanoTime();
			for(int i=0;i<l.size();i++) {
				 interation++;
//				 System.out.println(interation);
				 double time2 = System.nanoTime();
				 if((time2 - time1) / 1e9>10800) {
					 return lsol;
				 }
				 leaderloc=enumrate(fnum,P,l.get(i));
			     PSO pso=new PSO();
			     leaderattrac=pso.U_sol(instance, leaderloc, U, a,limit);
				 //System.out.println(Arrays.toString(leaderattrac));
				 FollowerSol fsol=subproblem(instance, leaderloc,leaderattrac,R,V,a,h, limit);
				 marketshare=calMarketShare(leaderattrac,fsol.getAttract(),instance,h);
				 if(marketshare>lsol.getLeaderMS()) {
					 lsol.setLeaderLoc(leaderloc);
					 lsol.setLeaderAttract(leaderattrac);
					 lsol.setFollowerLoc(fsol.getLocation());
					 lsol.setFollowerAttract(fsol.getAttract());
					 lsol.setLeaderMS(marketshare);
//					 System.out.println(Arrays.toString(leaderloc));
//					 System.out.println(Arrays.toString(leaderattrac));
//					 System.out.println(Arrays.toString(fsol.getLocation()));
//					 System.out.println(Arrays.toString(fsol.getAttract()));
//					 System.out.println(marketshare);
				}
			 }
			return lsol;
		}
		public static Solution solvePEH(Instance instance,double[] h,double[] a,int P,int R, double U, double V,double beta,double limit) throws IloException{
			Solution solution=new Solution();
			double time1 = System.nanoTime();
			int fnum=instance.getFnum();
			double[][] LeaderSol=new double[2][fnum];
			Sol sol=PEHAlg(instance,P,R,U,V,a,h, limit);		
			for(int i=0;i<fnum;i++) {
				LeaderSol[0][i]=sol.getLeaderLoc()[i];
			}
			LeaderSol[1]=sol.getLeaderAttract();
			
			double[][] FollwerSol=SOCP.SeparationSOCP(instance, LeaderSol, a, h, V, R, beta,limit);
			double time2 = System.nanoTime();                                   
			double Totaltime = (time2 - time1) / 1e9;//求解时间，单位s
			double obj=calMarketShare(LeaderSol[1],FollwerSol[1],instance,h);
			solution.setLeaderLoca(LeaderSol[0]);
			solution.setLeaderAttrac(LeaderSol[1]);
			solution.setFollowerLoca(FollwerSol[0]);
			solution.setFollowerAttrac(FollwerSol[1]);
			solution.setTime(Totaltime);
			solution.setObj(obj);
			solution.setInteration(interation);
			return solution;
			
		}
		 public static void main(String args[]) throws IOException, IloException{
		
			PEH p=new PEH();
			Main m=new Main();
			String pathname2="datasets/Facility_Leader.txt";
			File filename2 = new File(pathname2);
			String pathname3="datasets/Facility_Follower.txt";
			File filename3 = new File(pathname3);
			String pathname1="RandomData/"+"EX10/"+"100-10"+"/"+"100-10"+"-"+1+".txt";
			File filename1 = new File(pathname1);
			Instance instance=new Instance();
			instance=m.initData(filename1,filename2,filename3);
			double U=200;
			double V=200;
			int P=2;
			int R=2;
			int fnum=instance.getFnum();
			int cnum=instance.getCnum();
			double[] a=new double[fnum];
			for(int j=0;j<fnum;j++) {
				a[j]=50;
			}
			double[] h=new double[cnum];
			for(int i=0;i<cnum;i++) {
				h[i]=1.00/cnum;
			}
			double limit=0.1;
			double beta=1;
			double time1 = System.nanoTime();
			Solution sol=p.solvePEH(instance,h,a,P,R,U,V,beta,limit);
			double time2 = System.nanoTime();                                   
			double time = (time2 - time1) / 1e9;//求解时间，单位s
			
			System.out.println(Arrays.toString(sol.getLeaderLoca()));
			System.out.println(Arrays.toString(sol.getLeaderAttrac()));
			System.out.println(Arrays.toString(sol.getFollowerLoca()));
			System.out.println(Arrays.toString(sol.getFollowerAttrac()));
			System.out.println(sol.getObj());
			System.out.println(time);
			
//			List<List<Integer>> l=p.combine(fnum, P);
//			 for(int i=0;i<l.size();i++) {
//				 int[] x=enumrate(fnum,P,l.get(i));
//				 System.out.println(Arrays.toString(x));
//		 }
	}

}
