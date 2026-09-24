package Heuristic;

import java.util.Arrays;

public class Sol {
	//市场份额
    double LeaderMS;
    //选址
    int[] LeaderLoc;
    int[] FollowerLoc;
    //吸引力
    double[] LeaderAttract;
    double[] FollowerAttract;
    

    public Sol(double LeaderMS,int[] LeaderLoc,int[] FollowerLoc,double[] LeaderAttract,double[] FollowerAttract) {
		// TODO Auto-generated constructor stub
    	this.LeaderLoc=LeaderLoc;
    	this.FollowerLoc=FollowerLoc;
    	this.LeaderAttract=LeaderAttract;
    	this.FollowerAttract=FollowerAttract;
    	this.LeaderMS=LeaderMS;

	}

	public Sol copy() {
        return new Sol(LeaderMS, LeaderLoc.clone(),FollowerLoc.clone(),LeaderAttract.clone(),FollowerAttract.clone());
    }

    @Override
    public String toString() {
        return "Leader_marketshare = " + LeaderMS + " , Leader_location = " + Arrays.toString(LeaderLoc)+ " , Leader_attractiveness = " + Arrays.toString(LeaderAttract)+ " , Follower_location = " + Arrays.toString(FollowerLoc)+ " , Follower_attractiveness = " + Arrays.toString(FollowerAttract);
    }

	/**
	 * @return the leaderMS
	 */
	public double getLeaderMS() {
		return LeaderMS;
	}

	/**
	 * @param leaderMS the leaderMS to set
	 */
	public void setLeaderMS(double leaderMS) {
		LeaderMS = leaderMS;
	}

	/**
	 * @return the leaderLoc
	 */
	public int[] getLeaderLoc() {
		return LeaderLoc;
	}

	/**
	 * @param leaderLoc the leaderLoc to set
	 */
	public void setLeaderLoc(int[] leaderLoc) {
		LeaderLoc = leaderLoc;
	}

	/**
	 * @return the followerLoc
	 */
	public int[] getFollowerLoc() {
		return FollowerLoc;
	}

	/**
	 * @param followerLoc the followerLoc to set
	 */
	public void setFollowerLoc(int[] followerLoc) {
		FollowerLoc = followerLoc;
	}

	/**
	 * @return the leaderAttract
	 */
	public double[] getLeaderAttract() {
		return LeaderAttract;
	}

	/**
	 * @param leaderAttract the leaderAttract to set
	 */
	public void setLeaderAttract(double[] leaderAttract) {
		LeaderAttract = leaderAttract;
	}

	/**
	 * @return the followerAttract
	 */
	public double[] getFollowerAttract() {
		return FollowerAttract;
	}

	/**
	 * @param followerAttract the followerAttract to set
	 */
	public void setFollowerAttract(double[] followerAttract) {
		FollowerAttract = followerAttract;
	}
    
}
