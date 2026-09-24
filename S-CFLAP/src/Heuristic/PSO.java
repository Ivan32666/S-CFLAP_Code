package Heuristic;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import BC.Main;
import Common.Customer;
import Common.Facility;
import Common.Instance;

public class PSO {

    // 粒子对象
    class Particle {
        // 粒子速度数组（每个方向都有一个速度）
        double[] vArr;
        // 当前粒子坐标（自变量数组）
        double[] curVars;
        // 当前自变量对应的目标函数值
        double curObjValue;
        // 该粒子找到过的最佳目标函数值
        double bestObjValue;
        // 该粒子最好位置时的坐标
        double[] bestVars;
        // 全参构造
        public Particle(double[] vArr, double[] curVars, double curObjValue, double bestObjValue, double[] bestVars) {
            this.vArr = vArr;
            this.curVars = curVars;
            this.curObjValue = curObjValue;
            this.bestObjValue = bestObjValue;
            this.bestVars = bestVars;
        }
    }
    //惩罚因子
    double penalty=-1e3;
    // 粒子数量
    int n = 50;
    // 每个粒子的个体学习因子：自我认知，设置较大则不容易被群体带入局部最优，但会减缓收敛速度
    double c1 = 2;
    // 每个粒子的社会学习因子：社会认知，设置较大则加快收敛，但容易陷入局部最优
    double c2 = 2;
    // 粒子的惯性权重
    double w = 0.9;
    // 迭代的次数
    int MaxGen = 300;
    
    // 随机数对象
    Random random = new Random();
    
    // 粒子群
    Particle[] particles;
    // 最佳粒子
    Particle bestParticle;
    // 记录迭代过程
    public double[][][] positionArr;

    /**
     * @Description 初始化粒子的位置和速度
     */
    private void initParticles(double B, double[][] dis, ArrayList<Integer> list,int cnum,int varNum,double[] lbArr,double[] ubArr,double[] vMaxArr,double[] a,double limit) {
        // 初始化粒子群
        particles = new Particle[n];
        // 随机生成粒子
        for (int i = 0; i < particles.length; i++) {
            // 随机生成坐标和速度
            double[] vars = new double[varNum];
            double[] vArr = new double[varNum];
            for (int j = 0; j < varNum; j++) {
                vars[j] = random.nextDouble() * (ubArr[j] - lbArr[j]) + lbArr[j];
                vArr[j] = (random.nextDouble() - 0.5) * 2 * vMaxArr[j];
            }
            // 目标函数值
            double objValue = getObjValue(vars,B,dis,list,cnum,a,limit);
            particles[i] = new Particle(vArr.clone(), vars.clone(), objValue, objValue, vars.clone());
        }
        // 找到初始化粒子群中的最佳粒子
        bestParticle = copyParticle(particles[0]);
        for (int i = 1; i < particles.length; i++) {
            if (bestParticle.bestObjValue > particles[i].bestObjValue) {
                bestParticle = copyParticle(particles[i]);
            }
        }
    }

    /**
     * @Description 主要求解函数
     */
    public double[] solve(double B, double[][] dis, ArrayList<Integer> list,int cnum,int varNum,double[] lbArr,double[] ubArr,double[] vMaxArr,double[] a,double limit) {
        // 变量设置初步判断
        if(varNum != vMaxArr.length || varNum != lbArr.length || varNum != ubArr.length){
            throw new RuntimeException("变量维度不一致");
        }
        positionArr = new double[MaxGen][n][varNum];
        // 初始化粒子的位置和速度
        initParticles(B,dis,list,cnum,varNum,lbArr,ubArr,vMaxArr,a,limit);
        // 开始迭代
        for (int i = 0; i < MaxGen; i++) {
            // 依次更新第i个粒子的速度与位置
            for (int j = 0; j < particles.length; j++) {
                // 针对不同维度进行处理
                for (int k = 0; k < varNum; k++) {
                    // 更新速度
                    double newV = particles[j].vArr[k] * w
                            + c1 * random.nextDouble() * (particles[j].bestVars[k] - particles[j].curVars[k])
                            + c2 * random.nextDouble() * (bestParticle.bestVars[k] - particles[j].curVars[k]);
                    // 如果速度超过了最大限制，就对其进行调整
                    if (newV < -vMaxArr[k]) {
                        newV = -vMaxArr[k];
                    } else if (newV > vMaxArr[k]) {
                        newV = vMaxArr[k];
                    }
                    // 更新第j个粒子第k个维度上的位置
                    double newPos = particles[j].curVars[k] + newV;
                    // 记录迭代过程
                    positionArr[i][j][k] = newPos;
                    // 如果位置超出了定义域，就对其进行调整
                    if(newPos < lbArr[k]){
                        newPos = lbArr[k];
                    }else if(newPos > ubArr[k]){
                        newPos = ubArr[k];
                    }
                    // 赋值回去
                    particles[j].curVars[k] = newPos;
                    particles[j].vArr[k] = newV;
                }
                // 更新完所有维度后，再计算第j个粒子的函数值
                double objValueJ = getObjValue(particles[j].curVars,B,dis,list,cnum,a,limit);
                particles[j].curObjValue = objValueJ;
                if(objValueJ > particles[j].bestObjValue){
                    particles[j].bestVars = particles[j].curVars.clone();
                    particles[j].bestObjValue = particles[j].curObjValue;
                }
                if(objValueJ > bestParticle.bestObjValue){
                    bestParticle = copyParticle(particles[j]);
                }
            }
        }
        // 迭代结束，输出最优粒子位置和函数值
//        System.out.println("最优解为："+ bestParticle.bestObjValue);
//        System.out.println("最优解坐标为："+ Arrays.toString(bestParticle.bestVars));
        return bestParticle.bestVars;
    }

    /**
     * @param vars 自变量数组
     * @return 返回目标函数值
     */
    public double getObjValue(double[] vars,double B, double[][] dis, ArrayList<Integer> list,int cnum,double[] a,double limit) {
        
    	double p=0;
    	double sum=0;
    	double obj=0;
    	for(int i=0;i<vars.length;i++) {
    		sum=sum+vars[i];
    	}
    	p=p+Math.abs(sum-B);
    	
    	for(int j=0;j<vars.length;j++) {
    		if(vars[j]>0) {
    			p=p-Math.min(0, vars[j]-limit*B);
    		}
    	}
    	for(int i=0;i<cnum;i++) {
    		for(int j=0;j<list.size();j++) {
    			obj+=vars[j]/dis[i][list.get(j)];
    		}
    		
    	}
    	return penalty*p+obj;
        //return Math.pow(vars[0], 2) + Math.pow(vars[1], 2) - vars[0] * vars[1] - 10 * vars[0] - 4 * vars[1] + 60;
    }
    
    

    // 复制粒子
    public Particle copyParticle(Particle old) {
        return new Particle(old.vArr.clone(), old.curVars.clone(), old.curObjValue, old.bestObjValue, old.bestVars.clone());
    }
    
    public static double[][] dis(Instance instance){
		ArrayList<Customer> cus=instance.getCustomer();
		ArrayList<Facility> faC=instance.getFcandidate();
		int cnum=instance.getCnum();
		int fnum=instance.getFnum();
		double[][] d=new double[cnum][fnum];
		for(int i=0;i<cnum;i++) {
			for(int j=0;j<fnum;j++) {
				d[i][j]=(cus.get(i).getX()-faC.get(j).getX())*(cus.get(i).getX()-faC.get(j).getX())+(cus.get(i).getY()-faC.get(j).getY())*(cus.get(i).getY()-faC.get(j).getY()) ;	
			}
		}
		return d;
	}
    
    public double[] U_sol(Instance instance,int[] x,double B,double[] a,double limit) {
    	int cnum=instance.getCnum();
    	int fnum=instance.getFnum();
    	double[] u=new double[fnum];
    	ArrayList<Integer> list =new ArrayList<Integer>();
    	for(int j=0;j<x.length;j++) {
    		if(x[j]==1) {
    			list.add(j);
    		}
    	}
    	double[][] Dis=dis(instance);
    	int varNum = list.size();
	    // 自变量的上下界数组
	    double[] lbArr = new double[varNum];
	    double[] ubArr = new double[varNum];
	    // 粒子的每个维度上的最大速度（数组）
	    double[] vMaxArr = new double[varNum];
	    for(int i=0;i<varNum;i++) {
	    	lbArr[i]=0;
	    	ubArr[i]=a[i];
	    	vMaxArr[i]=1.2;
	    } 
		double[] u0=solve(B,Dis,list,cnum,varNum,lbArr,ubArr,vMaxArr,a,limit);
		int k=0;
		for(int j=0;j<fnum;j++) {
			if(x[j]==1) {
				u[j]=u0[k];
				k++;
			}
		}
    	return u;
    }
    
    public static void main(String args[]) throws IOException
    {
    	PSO pso=new PSO();
    	
    	Main m=new Main();
		String pathname2="datasets/Facility_Leader.txt";
		File filename2 = new File(pathname2);
		String pathname3="datasets/Facility_Follower.txt";
		File filename3 = new File(pathname3);
		String pathname1="RandomData/"+"EX1/"+"20-10"+"/"+"20-10"+"-"+1+".txt";
		File filename1 = new File(pathname1);
		Instance instance=new Instance();
		instance=m.initData(filename1,filename2,filename3);
		double B=100;
		double[][] Dis=dis(instance);
		int fnum=instance.getFnum();
		double[] a=new double[fnum];
		for(int j=0;j<fnum;j++) {
			a[j]=50;
		}
		double limit=0.2;
		int[] x= {0,1,1,0,0,1,1,0,0,0};
		double[] usol=pso.U_sol(instance,x,B,a,limit);
		System.out.println(Arrays.toString(usol));
    }
}

