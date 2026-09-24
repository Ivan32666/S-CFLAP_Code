package BC_Discrete;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;

import Common.Customer;
import Common.Facility;
import Common.Instance;
import Common.Solution;
import Heuristic.HeuristicAlg;
import Heuristic.PEH;
import ilog.concert.IloException;


public class Main {
		
	public static Instance initData(File filename1,File filename2,File filename3) throws IOException {
		Instance instance=new Instance();
		InputStreamReader read = new InputStreamReader(new FileInputStream(filename1));
        BufferedReader br = new BufferedReader(read);
        String lineTxt="" ;
        int line_location=0;
        int cnum=0;
    	int fnum=0;
    	ArrayList<Customer> customer=new ArrayList<Customer>();//顾客信息
    	ArrayList<Facility> Fcandidate=new ArrayList<Facility>();//候选设施点
    	ArrayList<Facility> FL=new ArrayList<Facility>();//leader已有设施
    	ArrayList<Facility> FF=new ArrayList<Facility>();//follower已有设施
    	
    	while ((lineTxt = br.readLine()) != null) {
        	line_location+=1;
        	if(line_location==1) {
        		cnum=Integer.parseInt(lineTxt.split(" ")[0]);
        		fnum=Integer.parseInt(lineTxt.split(" ")[1]);
        		instance.setCnum(cnum);
        		instance.setFnum(fnum);
        		
        	}
        	else if(line_location>1&&line_location<=cnum+1) {
        		Customer cus=new Customer();
        		cus.setX(Integer.parseInt(lineTxt.split(",")[0]));
        		cus.setY(Integer.parseInt(lineTxt.split(",")[1]));
        		customer.add(cus);
        	}else if(line_location>cnum+1) {
        		Facility fa=new Facility();    
        		fa.setX(Integer.parseInt(lineTxt.split(",")[0]));
        		fa.setY(Integer.parseInt(lineTxt.split(",")[1]));
        		Fcandidate.add(fa);
        	}
        }
    	instance.setCustomer(customer);
    	instance.setFcandidate(Fcandidate);

        InputStreamReader readFL = new InputStreamReader(new FileInputStream(filename2));
        BufferedReader brFL = new BufferedReader(readFL);
        String FLlineTxt="" ;
        int FLline_location=0;
        while ((FLlineTxt = brFL.readLine()) != null) {
        	FLline_location+=1;
        	if(FLline_location>1) {
        		Facility fa=new Facility();
        		fa.setX(Integer.parseInt(FLlineTxt.split(" ")[0]));
        		fa.setY(Integer.parseInt(FLlineTxt.split(" ")[1]));
        		fa.setSize(Double.parseDouble(FLlineTxt.split(" ")[2]));
        		FL.add(fa);
        	}
        }
        instance.setFL(FL);
        brFL.close();
        InputStreamReader readFF = new InputStreamReader(new FileInputStream(filename3));
        BufferedReader brFF = new BufferedReader(readFF);
        String FFlineTxt="" ;
        int FFline_location=0;
        while ((FFlineTxt = brFF.readLine()) != null) {
        	FFline_location+=1;
        	if(FFline_location>1) {
        		
        		Facility fa=new Facility();
        		fa.setX(Integer.parseInt(FFlineTxt.split(" ")[0]));
        		fa.setY(Integer.parseInt(FFlineTxt.split(" ")[1]));
        		fa.setSize(Double.parseDouble(FFlineTxt.split(" ")[2]));
        		FF.add(fa);
        	}
        }
        instance.setFF(FF);
        brFF.close();
        return instance;
	}
	
	
	public static String[] initDataName(File filename, int filenum) throws IOException {
		String[] name=new String[filenum];
		InputStreamReader read = new InputStreamReader(new FileInputStream(filename));
        BufferedReader br = new BufferedReader(read);
        int line_location=0;
        String line = br.readLine();
        while (line != null) {
        	name[line_location]=line;
        	line = br.readLine();
        	line_location+=1;
        }
        return name;
	}
	
	public static void batchSolve(int P,int R, double U, double V, int filenum, String namefilepath,double beta,double limit) throws IOException, IloException {
		
		
		File namefile = new File(namefilepath);
		String[] name=initDataName(namefile, filenum);
		String pathname2="datasets/Facility_Leader.txt";
		File filename2 = new File(pathname2);
		String pathname3="datasets/Facility_Follower.txt";
		File filename3 = new File(pathname3);
		
		
		for(int i=0;i<filenum;i++) {
//			ArrayList<Integer> list=new ArrayList<Integer>();
//
//			list.add(5);
//			list.add(1);
//
//
//			for(int j=0;j<list.size();j++) {
//				String pathname1="RandomData/"+"EX10/"+name[i]+"/"+name[i]+"-"+list.get(j)+".txt";
//				File filename1 = new File(pathname1);
//				Instance instance=new Instance();
//				instance=initData(filename1,filename2,filename3);
//				
//				int fnum=instance.getFnum();
//				int cnum=instance.getCnum();
//				double[] a= new double[fnum];
//				for(int k=0;k<fnum;k++) {
//					a[k]=50;
//				}
//				double[] h=new double[cnum];
//				for(int k=0;k<cnum;k++) {
//					h[k]=1.00/cnum;
//				}
//				
//				//Solution sol=staticCFLP.solveStaticCFLP(instance, h, a, P, R, U, V);
//				//Solution sol=UpAlg.solveOA(instance,h,a,P,R,U,V,beta,limit);
//				PEH peh=new PEH();
//				Solution sol=peh.solvePEH(instance,h,a,P,R,U,V,beta,limit);
//				CreateFile(fnum,cnum,P,R,U,V,sol,beta,list.get(j));
//			}
			ArrayList<Double> list1=new ArrayList<Double>();
			ArrayList<Double> list2=new ArrayList<Double>();
			ArrayList<Integer> list3=new ArrayList<Integer>();
			ArrayList<Double> list4=new ArrayList<Double>();
			
			for(int j=1;j<2;j++) {
				//String pathname1="RandomData/"+"EX9/"+name[i]+"/"+name[i]+"-"+j+".txt";
				String pathname1="RandomData/"+"EX for Table 7/"+name[i]+"-"+j+".txt";
				File filename1 = new File(pathname1);
				Instance instance=new Instance();
				instance=initData(filename1,filename2,filename3);
				
				
				int fnum=instance.getFnum();
				int cnum=instance.getCnum();
				double[] a= new double[fnum];
				for(int k=0;k<fnum;k++) {
					a[k]=100;
				}
				double[] h=new double[cnum];
				for(int k=0;k<cnum;k++) {
					h[k]=1.00/cnum;
				}
//				double U=a[0]*P/2; 
//				double V=a[0]*R/2;


				double[] u= {40,60,80,100};
				//Solution sol=staticCFLP.solveStaticCFLP(instance, h, a, P, R, U, V,limit);
				Solution sol=UpAlg_dis.solveOA(instance,h,a,u,P,R, U, V,beta,limit);
//				PEH peh=new PEH();
//				Solution sol=peh.solvePEH(instance,h,a,P,R,U,V,beta,limit);
				list1.add(sol.getTime());
				list2.add(sol.getSPtime());
				list3.add(sol.getCutnum());
				list4.add(sol.getObj());
//				list1.add(sol.getObj());
//				list2.add(sol.getTime());
				CreateFile(fnum,cnum,P,R,U,V,sol,beta,j);				
			}
			//RecordFile(name[i],list1,list2,list3,list4,P,R);
		}
		
	}
	
	public static void CreateFile(int fnum,int cnum,int P,int R,double U,double V,Solution solution,double beta,int j) throws IOException {

		//String fileName = "RandomResults/"+"EX9/"+cnum+"-"+fnum+"/"+cnum+"-"+fnum+"-"+P+"-"+R+"-"+U+"-"+V+"-"+j+"-"+"OA_Acc"+".txt";
		String fileName = "RandomResults/"+"EX_discrete/"+cnum+"-"+fnum+"-"+P+"-"+R+"-"+U+"-"+V+"-"+j+"-"+"Discrete1"+".txt";
		
		File file = new File(fileName);

		// 返回true表示文件成功

		// false 表示文件已经存在

		if (file.createNewFile()) {

			System.out.println("创建文件"+cnum+"-"+fnum+"-"+j+"成功！");

		} else {

			System.out.println("文件已经存在不需要重复创建");

		}

		// 使用FileWriter写文件

		try (FileWriter writer = new FileWriter(file)) {

			writer.write("Leader选址位置"+Arrays.toString(solution.getLeaderLoca()));
			writer.write("\r\n");
			writer.write("Leader选址规模"+Arrays.toString(solution.getLeaderAttrac()));
			writer.write("\r\n");
			writer.write("Follower选址位置"+Arrays.toString(solution.getFollowerLoca()));
			writer.write("\r\n");
			writer.write("Follower选址规模"+Arrays.toString(solution.getFollowerAttrac()));
			writer.write("\r\n"); 
			writer.write("运行时间"+solution.getTime());
			writer.write("\r\n");
			writer.write("Leader的市场占有率"+solution.getObj());
			writer.write("\r\n");
			writer.write("Cuts"+solution.getCutnum());
			writer.write("\r\n"); 
			writer.write("子问题运行时间"+solution.getSPtime());
			writer.write("\r\n"); 
//			writer.write("Gap:"+solution.getGap());
//			writer.write("\r\n");
			writer.write("迭代次数"+solution.getInteration());
		}
	}
	
	public static void RecordFile(String name,ArrayList<Double> list1,ArrayList<Double> list2,ArrayList<Integer> list3,ArrayList<Double> list4,int P,int R) throws IOException {

		String fileName = "RandomResults/"+"EX_discrete/"+"/"+name+"-"+P+"-"+R+".txt";
		//String fileName = "RandomResults/"+"Qi_OR/"+name+"/"+name+"-"+P+"-"+R+".txt";
		//String fileName = "RandomResults/"+cnum+"-"+fnum+"-"+P+"-"+R+"-"+U+"-"+V+"-"+j+".txt";
		
		File file = new File(fileName);

		// 返回true表示文件成功

		// false 表示文件已经存在

		if (file.createNewFile()) {

			System.out.println("创建文件成功！");

		} else {

			System.out.println("文件已经存在不需要重复创建");

		}

		// 使用FileWriter写文件

		try (FileWriter writer = new FileWriter(file)) {

			for(int i=0;i<list1.size();i++) {
				writer.write(""+list1.get(i));
				writer.write("\r\n");
			}
//			writer.write("子问题时间");
			writer.write("time");
			writer.write("\r\n");
			for(int i=0;i<list2.size();i++) {
				writer.write(""+list2.get(i));
				writer.write("\r\n");
			}
			writer.write("cut");
			writer.write("\r\n");
			for(int i=0;i<list3.size();i++) {
				writer.write(""+list3.get(i));
				writer.write("\r\n");
			}
			writer.write("目标函数值");
			writer.write("\r\n");
			for(int i=0;i<list4.size();i++) {
				writer.write(""+list4.get(i));
				writer.write("\r\n");
			}
		}
	}

	public static void main(String[] args) throws IOException, IloException { 
		int P=3;
		int R=3;
		int filenum=1;
		String namefilepath1="RandomData/data.txt";
		double beta=2;
		double limit=0.4;
		double U=150; 
		double V=150;
		batchSolve(P,R,U,V,filenum,namefilepath1,beta,limit);
		
//		filenum=1;
//		P=3;R=3;
//		U=150;V=150;
//		batchSolve(P,R,U,V,filenum,namefilepath1,beta,limit);
//		P=3;R=3;
//		batchSolve(P,R,filenum,namefilepath1,beta,limit);
//		
	}
}
