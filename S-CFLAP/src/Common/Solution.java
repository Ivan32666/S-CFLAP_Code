package Common;

public class Solution {
	private double[] leaderLoca;
	private double[] leaderAttrac;
	private double[] followerLoca;
	private double[] followerAttrac;
	private double time;
	private double obj;
	private int cutnum;
	private double SPtime;
	private double Gap;
	private double Interation;
	/**
	 * @return the interation
	 */
	public double getInteration() {
		return Interation;
	}
	/**
	 * @param interation the interation to set
	 */
	public void setInteration(double interation) {
		Interation = interation;
	}
	/**
	 * @return the gap
	 */
	public double getGap() {
		return Gap;
	}
	/**
	 * @param gap the gap to set
	 */
	public void setGap(double gap) {
		Gap = gap;
	}
	/**
	 * @return the sPtime
	 */
	public double getSPtime() {
		return SPtime;
	}
	/**
	 * @param sPtime the sPtime to set
	 */
	public void setSPtime(double sPtime) {
		SPtime = sPtime;
	}
	/**
	 * @return the cutnum
	 */
	public int getCutnum() {
		return cutnum;
	}
	/**
	 * @param cutnum the cutnum to set
	 */
	public void setCutnum(int cutnum) {
		this.cutnum = cutnum;
	}
	/**
	 * @return the obj
	 */
	public double getObj() {
		return obj;
	}
	/**
	 * @param obj the obj to set
	 */
	public void setObj(double obj) {
		this.obj = obj;
	}
	/**
	 * @return the leaderLoca
	 */
	public double[] getLeaderLoca() {
		return leaderLoca;
	}
	/**
	 * @param leaderLoca the leaderLoca to set
	 */
	public void setLeaderLoca(double[] leaderLoca) {
		this.leaderLoca = leaderLoca;
	}
	/**
	 * @return the leaderAttrac
	 */
	public double[] getLeaderAttrac() {
		return leaderAttrac;
	}
	/**
	 * @param leaderAttrac the leaderAttrac to set
	 */
	public void setLeaderAttrac(double[] leaderAttrac) {
		this.leaderAttrac = leaderAttrac;
	}
	/**
	 * @return the followerLoca
	 */
	public double[] getFollowerLoca() {
		return followerLoca;
	}
	/**
	 * @param followerLoca the followerLoca to set
	 */
	public void setFollowerLoca(double[] followerLoca) {
		this.followerLoca = followerLoca;
	}
	/**
	 * @return the followerAttrac
	 */
	public double[] getFollowerAttrac() {
		return followerAttrac;
	}
	/**
	 * @param followerAttrac the followerAttrac to set
	 */
	public void setFollowerAttrac(double[] followerAttrac) {
		this.followerAttrac = followerAttrac;
	}
	/**
	 * @return the time
	 */
	public double getTime() {
		return time;
	}
	/**
	 * @param time the time to set
	 */
	public void setTime(double time) {
		this.time = time;
	}
	
}
