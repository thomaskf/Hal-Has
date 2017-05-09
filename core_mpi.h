/*
 *
 * core_mpi.h
 * HAL_HAS
 *
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2014, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * All rights reserved. CSIRO is willing to grant you a license to HAL-HAS on the terms of the GNU General Public
 * License version 3 as published by the Free Software Foundation (http://www.gnu.org/licenses/gpl.html), except
 * where otherwise indicated for third party material.
 * The following additional terms apply under clause 7 of that license:
 * EXCEPT AS EXPRESSLY STATED IN THIS AGREEMENT AND TO THE FULL EXTENT PERMITTED BY APPLICABLE LAW, THE SOFTWARE
 * IS PROVIDED "AS-IS". CSIRO MAKES NO REPRESENTATIONS, WARRANTIES OR CONDITIONS OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO ANY REPRESENTATIONS, WARRANTIES OR CONDITIONS REGARDING THE CONTENTS OR ACCURACY
 * OF THE SOFTWARE, OR OF TITLE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, THE ABSENCE
 * OF LATENT OR OTHER DEFECTS, OR THE PRESENCE OR ABSENCE OF ERRORS, WHETHER OR NOT DISCOVERABLE.
 * TO THE FULL EXTENT PERMITTED BY APPLICABLE LAW, IN NO EVENT SHALL CSIRO BE LIABLE ON ANY LEGAL THEORY (INCLUDING,
 * WITHOUT LIMITATION, IN AN ACTION FOR BREACH OF CONTRACT, NEGLIGENCE OR OTHERWISE) FOR ANY CLAIM, LOSS, DAMAGES
 * OR OTHER LIABILITY HOWSOEVER INCURRED.  WITHOUT LIMITING THE SCOPE OF THE PREVIOUS SENTENCE THE EXCLUSION OF
 * LIABILITY SHALL INCLUDE: LOSS OF PRODUCTION OR OPERATION TIME, LOSS, DAMAGE OR CORRUPTION OF DATA OR RECORDS;
 * OR LOSS OF ANTICIPATED SAVINGS, OPPORTUNITY, REVENUE, PROFIT OR GOODWILL, OR OTHER ECONOMIC LOSS; OR ANY SPECIAL,
 * INCIDENTAL, INDIRECT, CONSEQUENTIAL, PUNITIVE OR EXEMPLARY DAMAGES, ARISING OUT OF OR IN CONNECTION WITH THIS
 * AGREEMENT, ACCESS OF THE SOFTWARE OR ANY OTHER DEALINGS WITH THE SOFTWARE, EVEN IF CSIRO HAS BEEN ADVISED OF
 * THE POSSIBILITY OF SUCH CLAIM, LOSS, DAMAGES OR OTHER LIABILITY.
 * APPLICABLE LEGISLATION SUCH AS THE AUSTRALIAN CONSUMER LAW MAY APPLY REPRESENTATIONS, WARRANTIES, OR CONDITIONS,
 * OR IMPOSES OBLIGATIONS OR LIABILITY ON CSIRO THAT CANNOT BE EXCLUDED, RESTRICTED OR MODIFIED TO THE FULL EXTENT
 * SET OUT IN THE EXPRESS TERMS OF THIS CLAUSE ABOVE "CONSUMER GUARANTEES".  TO THE EXTENT THAT SUCH CONSUMER
 * GUARANTEES CONTINUE TO APPLY, THEN TO THE FULL EXTENT PERMITTED BY THE APPLICABLE LEGISLATION, THE LIABILITY
 * OF CSIRO UNDER THE RELEVANT CONSUMER GUARANTEE IS LIMITED (WHERE PERMITTED AT CSIRO’S OPTION) TO ONE OF FOLLOWING
 * REMEDIES OR SUBSTANTIALLY EQUIVALENT REMEDIES:
 * (a)               THE REPLACEMENT OF THE SOFTWARE, THE SUPPLY OF EQUIVALENT SOFTWARE, OR SUPPLYING RELEVANT
 *                   SERVICES AGAIN;
 * (b)               THE REPAIR OF THE SOFTWARE;
 * (c)               THE PAYMENT OF THE COST OF REPLACING THE SOFTWARE, OF ACQUIRING EQUIVALENT SOFTWARE, HAVING THE
 *                   RELEVANT SERVICES SUPPLIED AGAIN, OR HAVING THE SOFTWARE REPAIRED.
 * IN THIS CLAUSE, CSIRO INCLUDES ANY THIRD PARTY AUTHOR OR OWNER OF ANY PART OF THE SOFTWARE OR MATERIAL DISTRIBUTED
 * WITH IT.  CSIRO MAY ENFORCE ANY RIGHTS ON BEHALF OF THE RELEVANT THIRD PARTY.
 * Third Party Components
 * The following third party components are distributed with the Software.  You agree to comply with the license
 * terms for these components as part of accessing the Software.  Other third party software may also be identified
 * in separate files distributed with the Software.
 * ___________________________________________________________________
 * 
 * R : A Computer Language for Statistical Data Analysis version 3.0.1 (http://cran.r-project.org/src/base/R-3/R-3.0.1.tar.gz)
 * Copyright (C) 2000-2004 The R Core Team
 * This software is licensed under GNU GPL
 * 
 * JACOBI_EIGENVALUE.C (http://people.sc.fsu.edu/~jburkardt/c_src/jacobi_eigenvalue/jacobi_eigenvalue.c)
 * Copyright (C) 2003-2013 John Burkardt
 * This software is licensed under GNU LGPL (http://www.gnu.org/licenses/lgpl.html)
 * ___________________________________________________________________
 */


#ifndef __RAL_RAS__CORE_MPI__
#define __RAL_RAS__CORE_MPI__

#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include "parameters.h"
#include "core.h"

using namespace std;

// ==================================================================================================
// functions for MPI
// ==================================================================================================

// for HAL info
void sendInfo(vector<int>* rateMatrix , int jobID, int currMode, int currMaxIT, int proceed,
		int num_edges, int num_w, int num_chars, int destID, int tagID, ParameterSet* ps, VariableSet* vs);

void receiveInfo(vector<int>& rateMatrix , int& jobID, int& currMode, int& currMaxIT, int& num_rateMat, int& proceed, 
		int num_edges, int num_w, int num_chars, int tagID, ParameterSet* ps, VariableSet* vs, int& isPsVsUpdated);

void sendPsVs(ParameterSet& ps, VariableSet& vs, int destID, int tagID);

void receivePsVs(ParameterSet& ps, VariableSet& vs, int tagID);

/*
void sendInfo(vector<int>* rateMatrix , int jobID, int currMode, int currMaxIT, int num_rateMat, int signal,
		int num_edges, int destID, int tagID);
		
void receiveInfo(vector<int>& rateMatrix , int& jobID, int& currMode, int& currMaxIT, int& num_rateMat, int& signal,
		int num_edges, int tagID);

*/

// for HAL result
void sendResult(int jobID, double logL, double ic, ParameterSet& ps, VariableSet& vs, int numIter, int tagID);

void receiveResult(int& jobID, double& logL, double& ic, ParameterSet& ps, VariableSet& vs,
					int& numIter, int tagID, int& sender);

// for HAS info
void sendInfoHAS(vector<int>* rateMatrix , int jobID, int currMode, int currMaxIT, int num_rateMat, int proceed,
		int numRateCat, int modelType, int num_edges, int destID, int tagID);

void receiveInfoHAS(vector<int>& rateMatrix , int& jobID, int& currMode, int& currMaxIT, int& num_rateMat, int& proceed,
		int& numRateCat, int& modelType, int num_edges, int tagID);

// for HAS result
void sendResultHAS(int jobID, double logL, double ic, int df, AllParameterSet& all_ps, VariableSet& vs, int numIter, int tagID);

void receiveResultHAS(int& jobID, double& logL, double& ic, int& df, AllParameterSet& all_ps, VariableSet& vs,
					int& numIter, int tagID, int& sender);

// Terminate all processes
void terminateAllProcs(int num_edges, int num_w, int num_chars, int tagID);

void terminateAllProcsHAS(int num_edges, int tagID);

#endif /* defined(__RAL_RAS__CORE_MPI__) */