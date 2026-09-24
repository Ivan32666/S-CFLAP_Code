package Common;

import java.util.ArrayList;

import ilog.concert.IloNumExpr;

public class Cuts {
	ArrayList<IloNumExpr> lhs;
	ArrayList<Double> rhs;
	/**
	 * @return the lhs
	 */
	public ArrayList<IloNumExpr> getLhs() {
		return lhs;
	}
	/**
	 * @param lhs the lhs to set
	 */
	public void setLhs(ArrayList<IloNumExpr> lhs) {
		this.lhs = lhs;
	}
	/**
	 * @return the rhs
	 */
	public ArrayList<Double> getRhs() {
		return rhs;
	}
	/**
	 * @param rhs the rhs to set
	 */
	public void setRhs(ArrayList<Double> rhs) {
		this.rhs = rhs;
	}

}
