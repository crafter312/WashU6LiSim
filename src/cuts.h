/* Header file containing 2D cuts for analysis purposes
 * Written 9 April 2026 by Henry Webb (h.s.webb@wustl.edu)
 * 
 * Static initialization method by ROOT forum user Axel, May 2015
 * https://root-forum.cern.ch/t/root-6-02-08-loading-and-using-cuts/19219/6
 * 
 * I don't know if SetVarX/Y are needed when using this in code,
 * put I included them anyways. At the very least, this allows
 * you to look at the cuts and understand how to use them (or
 * at least what variables they should be used with).
 */

#include <TCutG.h>

/******** CUTS FOR REACTION DEPTH RECONSTRUCTION ********/

using double_array = double[];

// 13C ground state
TCutG cut_gs("cut_gs", 7, (double_array){41.796475, 28.058288, 30.77167860421347, 45.84021747609521, 52.90030212236147, 50.8981885659576, 41.796475}, (double_array){8.6033205, 22.350447, 23.37982134326133, 9.79872280116815, 5.481992213216774, 4.087356177109406, 8.6033205});
struct cut_gs_Init_t {
	cut_gs_Init_t() {
		cut_gs.SetVarX("reconFragments[0].reconEnergy+reconFragments[1].reconEnergy");
		cut_gs.SetVarY("targEloss");
	}
} cut_gs_Init;

// 13C 5/2+ state
TCutG cut_52p("cut_52p", 8, (double_array){41.796475, 28.058288, 25.79273831263017, 24.55458913959094, 37.91079404875884, 46.49880746175437, 50.8981885659576, 41.796475}, (double_array){8.6033205, 22.350447, 23.31341010344669, 23.31341010344669, 7.208684448397324, 4.917496674792361, 4.087356177109406, 8.6033205});
struct cut_52p_Init_t {
	cut_52p_Init_t() {
		cut_52p.SetVarX("reconFragments[0].reconEnergy+reconFragments[1].reconEnergy");
		cut_52p.SetVarY("targEloss");
	}
} cut_52p_Init;
