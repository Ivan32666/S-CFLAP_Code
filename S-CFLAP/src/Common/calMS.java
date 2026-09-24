package Common;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import com.mathworks.toolbox.javabuilder.MWException;

import BC.Main;
import BC.OA;
import BC.UpAlg;
import ilog.concert.IloException;

public class calMS {
	
	public static double calDistanceSquare(Customer cus,Facility fa) {
		double ds=0;
		ds=(cus.getX()-fa.getX())*(cus.getX()-fa.getX())+(cus.getY()-fa.getY())*(cus.getY()-fa.getY());	
		return ds;
	}
	
	
	public static double calMarketShare(double[] x,double[] y,Instance instance,double[] u,double[] v,double[] h,double beta) {
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
			share+= h[i]*shareL/(shareL+shareF);
			
		}

		return share;
	}
	
public static void main(String[] args) throws IOException, IloException, MWException{
		
		Main m=new Main();
		String pathname2="datasets/Facility_Leader.txt";
		File filename2 = new File(pathname2);
		String pathname3="datasets/Facility_Follower.txt";
		File filename3 = new File(pathname3);
		String pathname1="RandomData/"+"Qi_OR/"+"100-40"+".txt";
		File filename1 = new File(pathname1);
		Instance instance=new Instance();
		instance=m.initData(filename1,filename2,filename3);
	
		int cnum=instance.getCnum();
		int fnum=instance.getCnum();
		int R=2;
		int P=2;
		double V=80;
		double U=80;
		double[] h=new double[cnum];
		for(int k=0;k<cnum;k++) {
			h[k]=1.00/cnum;
		}
		double[] a= new double[fnum];
		for(int j=0;j<fnum;j++) {
			a[j]=50;
		}
		double limit=0.4;
		double beta=2;
		double[][] LeaderSol= {
				{0.0, -0.0, 0.0, 0.0, 0.0, 0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0,0.0, -0.0, 0.0, 0.0, 0.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 1.0, -0.0, 0.0},{0.0, -0.0, 0.0, 0.0, 0.0, 0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0,0.0, -0.0, 0.0, 0.0, 0.0, 40.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 40.0, -0.0, 0.0}
};
		double[][] FollwerSol=OA.SeparationProblem(LeaderSol, instance, V, a, h, R,limit,beta);
		System.out.println(Arrays.toString(FollwerSol[0]));
		System.out.println(Arrays.toString(FollwerSol[1]));
		double obj=calMarketShare(LeaderSol[0],FollwerSol[0],instance,LeaderSol[1],FollwerSol[1],h,beta);


		System.out.println(obj);
		
		
	}
}
