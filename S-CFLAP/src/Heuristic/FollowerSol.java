package Heuristic;

import java.util.Arrays;

public class FollowerSol {
	//市场份额
    double marketshare;
    //选址
    int[] location;
    //吸引力
    double[] attract;
    

    public FollowerSol(double marketshare, int[] location,double[] attract) {
		// TODO Auto-generated constructor stub
    	this.location=location;
    	this.marketshare=marketshare;
    	this.attract=attract;
	}

	public FollowerSol copy() {
        return new FollowerSol(marketshare, location.clone(),attract);
    }

    @Override
    public String toString() {
        return "marketshare = " + marketshare + " , location = " + Arrays.toString(location)+ " , design = " + Arrays.toString(attract);
    }

	/**
	 * @return the attract
	 */
	public double[] getAttract() {
		return attract;
	}

	/**
	 * @param attract the attract to set
	 */
	public void setAttract(double[] attract) {
		this.attract = attract;
	}

	/**
	 * @return the marketshare
	 */
	public double getMarketshare() {
		return marketshare;
	}

	/**
	 * @param marketshare the marketshare to set
	 */
	public void setMarketshare(double marketshare) {
		this.marketshare = marketshare;
	}

	/**
	 * @return the location
	 */
	public int[] getLocation() {
		return location;
	}

	/**
	 * @param location the location to set
	 */
	public void setLocation(int[] location) {
		this.location = location;
	}

}
