package Common;

import java.util.ArrayList;

import Common.Customer;
import Common.Facility;

public class Instance {
	private ArrayList<Customer> customer=new ArrayList<Customer>();//顾客信息
	private ArrayList<Facility> Fcandidate=new ArrayList<Facility>();//候选设施点
	private ArrayList<Facility> FL=new ArrayList<Facility>();//leader已有设施
	private ArrayList<Facility> FF=new ArrayList<Facility>();//follower已有设施
	private int fnum;
	private int cnum;
	/**
	 * @return the customer
	 */
	public ArrayList<Customer> getCustomer() {
		return customer;
	}
	/**
	 * @param customer the customer to set
	 */
	public void setCustomer(ArrayList<Customer> customer) {
		this.customer = customer;
	}
	/**
	 * @return the fcandidate
	 */
	public ArrayList<Facility> getFcandidate() {
		return Fcandidate;
	}
	/**
	 * @param fcandidate the fcandidate to set
	 */
	public void setFcandidate(ArrayList<Facility> fcandidate) {
		Fcandidate = fcandidate;
	}
	/**
	 * @return the fL
	 */
	public ArrayList<Facility> getFL() {
		return FL;
	}
	/**
	 * @param fL the fL to set
	 */
	public void setFL(ArrayList<Facility> fL) {
		FL = fL;
	}
	/**
	 * @return the fF
	 */
	public ArrayList<Facility> getFF() {
		return FF;
	}
	/**
	 * @param fF the fF to set
	 */
	public void setFF(ArrayList<Facility> fF) {
		FF = fF;
	}
	/**
	 * @return the fnum
	 */
	public int getFnum() {
		return fnum;
	}
	/**
	 * @param fnum the fnum to set
	 */
	public void setFnum(int fnum) {
		this.fnum = fnum;
	}
	/**
	 * @return the cnum
	 */
	public int getCnum() {
		return cnum;
	}
	/**
	 * @param cnum the cnum to set
	 */
	public void setCnum(int cnum) {
		this.cnum = cnum;
	}
	
	
	
}
