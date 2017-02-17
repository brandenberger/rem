//####################################################################include <string.h>
#include <Rcpp.h>
using namespace Rcpp;
  
//TODO: tidy up functions - within 80char/line
//TODO: similarity-Average-Total: why does target-sim not have "match"-option for vector w?
  
//####################################################################
//####################################################################
//####################################################################

//####################################################################
// [[Rcpp::export]]
NumericVector inertiaCpp(
  NumericVector time,
  NumericVector weightvar,
  std::vector<std::string> sender,
  std::vector<std::string> target,
  std::vector<std::string> typevar,
  std::string type1, 
  std::string type2, 
  std::vector<std::string> attrvar,
  std::string attr1,
  std::string attr2, 
  double xlog, 
  std::string inertiatype ) {
  
  NumericVector result(sender.size());
  
  // for-loop i: for each event, do:
    for ( int i = 0; i < sender.size(); i++){
      
      // reset all the variables
      double inertia = 0;
      double totalinertia = 0;
      double weight = 0;
      double resulttemp = 0;
      
      // for-loop w: go back over all events and filter them
      for ( int w = 0; w < i; w++ ){
        
        // sender-target inertia only
        if ( inertiatype == "s-t-only" ){
          if ( target[w] == target[i] && sender[w] == sender[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typematch 
        if ( inertiatype == "s-t-typematch" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[w] == typevar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typemix
        if ( inertiatype == "s-t-typemix" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[i] == type1 && typevar[w] == type2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typefilter
        if ( inertiatype == "s-t-typefilter" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[w] == type1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-attributematch
        if ( inertiatype == "s-t-attributematch" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-attributemix
        if ( inertiatype == "s-t-attributemix" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-attributefilter
        if ( inertiatype == "s-t-attributefilter" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typematch-attributematch
        if ( inertiatype == "s-t-typematch-attributematch" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[w] == typevar[i] && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typematch-attributemix
        if ( inertiatype == "s-t-typematch-attributemix" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[w] == typevar[i] && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typematch-attributefilter
        if ( inertiatype == "s-t-typematch-attributefilter" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[w] == typevar[i] && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typemix-attributematch
        if ( inertiatype == "s-t-typemix-attributematch" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typemix-attributemix
        if ( inertiatype == "s-t-typemix-attributemix" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typemix-attributefilter
        if ( inertiatype == "s-t-typemix-attributefilter" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typefilter-attributematch
        if ( inertiatype == "s-t-typefilter-attributematch" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[w] == type1 && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typefilter-attributemix
        if ( inertiatype == "s-t-typefilter-attributemix" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[w] == type1 && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typefilter-attributefilter
        if ( inertiatype == "s-t-typefilter-attributefilter" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && typevar[w] == type1 && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        //
          // calculate the weight
        if ( time[i] == time[w] ){
          totalinertia = totalinertia + 0;
        } else {
          inertia = weight * exp( - ( time[i] - time[w] ) * xlog) * xlog ;
          totalinertia = totalinertia + inertia;
        }
        
        resulttemp = totalinertia;
      } // closes w-loop
      
      // ascribe calculated weights to the result-variable		
      result[i] = resulttemp;
      
    } // closes i-loop
  
  // return variable and hand it back to R
  return Rcpp::wrap(result);
}


//####################################################################
// [[Rcpp::export]]
NumericVector degreeCpp(
  NumericVector time,
  NumericVector weightvar,
  std::vector<std::string> degreevar,
  std::vector<std::string> typevar,
  std::string type1, 
  std::string type2, 
  std::vector<std::string> attrvar,
  std::string attr1,
  std::string attr2, 
  double xlog, 
  std::string degreetype ) {
  
  NumericVector result(degreevar.size());
  
  // for-loop i: for each event, do:
    for ( int i = 0; i < degreevar.size(); i++){
      
      // reset all the variables
      double degree = 0;
      double totaldegree = 0;
      double weight = 0;
      double resulttemp = 0;
      
      // for-loop w: go back over all events and filter them
      for ( int w = 0; w < i; w++ ){
        
        // degreevar only
        if ( degreetype == "d-only" ){
          if ( degreevar[w] == degreevar[i]  ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typematch 
        if ( degreetype == "d-typematch" ){
          if ( degreevar[w] == degreevar[i] && typevar[w] == typevar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typemix
        if ( degreetype == "d-typemix" ){
          if ( degreevar[w] == degreevar[i] && typevar[i] == type1 && typevar[w] == type2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typefilter
        if ( degreetype == "d-typefilter" ){
          if ( degreevar[w] == degreevar[i] && typevar[w] == type1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-attributematch
        if ( degreetype == "d-attributematch" ){
          if ( degreevar[w] == degreevar[i] && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-attributemix
        if ( degreetype == "d-attributemix" ){
          if ( degreevar[w] == degreevar[i] && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-attributefilter
        if ( degreetype == "d-attributefilter" ){
          if ( degreevar[w] == degreevar[i] && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typematch-attributematch
        if ( degreetype == "d-typematch-attributematch" ){
          if ( degreevar[w] == degreevar[i] && typevar[w] == typevar[i] && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typematch-attributemix
        if ( degreetype == "d-typematch-attributemix" ){
          if ( degreevar[w] == degreevar[i] && typevar[w] == typevar[i] && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typematch-attributefilter
        if ( degreetype == "d-typematch-attributefilter" ){
          if ( degreevar[w] == degreevar[i] && typevar[w] == typevar[i] && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typemix-attributematch
        if ( degreetype == "d-typemix-attributematch" ){
          if ( degreevar[w] == degreevar[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typemix-attributemix
        if ( degreetype == "d-typemix-attributemix" ){
          if ( degreevar[w] == degreevar[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typemix-attributefilter
        if ( degreetype == "d-typemix-attributefilter" ){
          if ( degreevar[w] == degreevar[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typefilter-attributematch
        if ( degreetype == "d-typefilter-attributematch" ){
          if ( degreevar[w] == degreevar[i] && typevar[w] == type1 && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typefilter-attributemix
        if ( degreetype == "d-typefilter-attributemix" ){
          if ( degreevar[w] == degreevar[i] && typevar[w] == type1 && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typefilter-attributefilter
        if ( degreetype == "d-typefilter-attributefilter" ){
          if ( degreevar[w] == degreevar[i] && typevar[w] == type1 && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        //
          // calculate the weight
        if ( time[i] == time[w] ){
          totaldegree = totaldegree + 0;
        } else {
          degree = weight * exp( - ( time[i] - time[w] ) * xlog) * xlog ;
          totaldegree = totaldegree + degree;
        }
        
        resulttemp = totaldegree;
      } // closes w-loop
      
      // ascribe calculated weights to the result-variable		
      result[i] = resulttemp;
      
    } // closes i-loop
  
  // return variable and hand it back to R
  return Rcpp::wrap(result);
}

//####################################################################
// [[Rcpp::export]]
NumericVector degreeOneModeCpp(
  NumericVector time,
  NumericVector weightvar,
  std::vector<std::string> degreevar,
  std::vector<std::string> degreeothermodevar,
  std::vector<std::string> typevar,
  std::string type1, 
  std::string type2, 
  std::vector<std::string> attrvar,
  std::string attr1,
  std::string attr2, 
  double xlog, 
  std::string degreetype ) {
  
  NumericVector result(degreevar.size());
  
  // for-loop i: for each event, do:
    for ( int i = 0; i < degreevar.size(); i++){
      
      // reset all the variables
      double degree = 0;
      double totaldegree = 0;
      double weight = 0;
      double resulttemp = 0;
      
      // for-loop w: go back over all events and filter them
      for ( int w = 0; w < i; w++ ){
        
        // degreevar only
        if ( degreetype == "d-only" ){
          if ( degreeothermodevar[w] == degreevar[i]  ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typematch 
        if ( degreetype == "d-typematch" ){
          if ( degreeothermodevar[w] == degreevar[i] && typevar[w] == typevar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typemix
        if ( degreetype == "d-typemix" ){
          if ( degreeothermodevar[w] == degreevar[i] && typevar[i] == type1 && typevar[w] == type2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typefilter
        if ( degreetype == "d-typefilter" ){
          if ( degreeothermodevar[w] == degreevar[i] && typevar[w] == type1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-attributematch
        if ( degreetype == "d-attributematch" ){
          if ( degreeothermodevar[w] == degreevar[i] && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-attributemix
        if ( degreetype == "d-attributemix" ){
          if ( degreeothermodevar[w] == degreevar[i] && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-attributefilter
        if ( degreetype == "d-attributefilter" ){
          if ( degreeothermodevar[w] == degreevar[i] && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typematch-attributematch
        if ( degreetype == "d-typematch-attributematch" ){
          if ( degreeothermodevar[w] == degreevar[i] && typevar[w] == typevar[i] && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typematch-attributemix
        if ( degreetype == "d-typematch-attributemix" ){
          if ( degreeothermodevar[w] == degreevar[i] && typevar[w] == typevar[i] && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typematch-attributefilter
        if ( degreetype == "d-typematch-attributefilter" ){
          if ( degreeothermodevar[w] == degreevar[i] && typevar[w] == typevar[i] && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typemix-attributematch
        if ( degreetype == "d-typemix-attributematch" ){
          if ( degreeothermodevar[w] == degreevar[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typemix-attributemix
        if ( degreetype == "d-typemix-attributemix" ){
          if ( degreeothermodevar[w] == degreevar[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typemix-attributefilter
        if ( degreetype == "d-typemix-attributefilter" ){
          if ( degreeothermodevar[w] == degreevar[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typefilter-attributematch
        if ( degreetype == "d-typefilter-attributematch" ){
          if ( degreeothermodevar[w] == degreevar[i] && typevar[w] == type1 && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typefilter-attributemix
        if ( degreetype == "d-typefilter-attributemix" ){
          if ( degreeothermodevar[w] == degreevar[i] && typevar[w] == type1 && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // degree-typefilter-attributefilter
        if ( degreetype == "d-typefilter-attributefilter" ){
          if ( degreeothermodevar[w] == degreevar[i] && typevar[w] == type1 && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        //
          // calculate the weight
        if ( time[i] == time[w] ){
          totaldegree = totaldegree + 0;
        } else {
          degree = weight * exp( - ( time[i] - time[w] ) * xlog) * xlog ;
          totaldegree = totaldegree + degree;
        }
        
        resulttemp = totaldegree;
      } // closes w-loop
      
      // ascribe calculated weights to the result-variable		
      result[i] = resulttemp;
      
    } // closes i-loop
  
  // return variable and hand it back to R
  return Rcpp::wrap(result);
}

//####################################################################
// [[Rcpp::export]]
double fourCycleCpp(
  std::vector<std::string> sender,
  std::string currentSender,
  std::vector<std::string> target,
  std::string currentTarget,
  std::vector<std::string> typevar,
  std::string currentType,
  NumericVector time,
  double currentTime,
  NumericVector weightvar,
  double xlog,
  std::vector<std::string> attrvarAaj,
  std::string attrAaj, 
  std::vector<std::string> attrvarBib,
  std::string attrBib,
  std::vector<std::string> attrvarCij,
  std::string attrCij, 
  std::string fourCycleType, 
  std::vector<std::string> w, //what else has a said?
  std::vector<std::string> x, //who else has used a (same opinion = positive, opposite oppingion = negative)
  int i, 
  int begin) {
  
  double result;
  std::vector<std::string> y;
  std::vector<std::string> wy;
  double weightA;
  double weightB;
  double weightC;
  double tempTotalWeightC;
  double tempTotalWeightA;
  double tempTotalWeightAPositive;
  double tempTotalWeightANegative;
  double tempTotalWeightB;
  double tempTotalWeightCPositive;
  double tempTotalWeightCNegative;
  double tempTotalWeightABC;
  double tempTotalWeightABCPositive;
  double tempTotalWeightABCNegative;
  double totalWeightABCPositive;
  double totalWeightABCNegative;
  double totalWeightABC;
  
  //in R-loop now: Filter: only current events with given attribute are selected.
  //if ( attrvarNow[i] == attrNow ) {
    
    y.clear();
    totalWeightABCPositive = 0;
    totalWeightABCNegative = 0;
    tempTotalWeightABC = 0;
    tempTotalWeightABCNegative = 0;
    tempTotalWeightABCPositive = 0;
    totalWeightABC = 0;
    
    // for each person in the list x (m-loop open) (=y-vector; wy-vector)
      for (int m = 0; m < x.size(); m++ ) {
        
        tempTotalWeightB = 0;
        
        // What did actor say in past? (n-loop open)
        y.clear();
        for ( int n = begin; n < i; n++ ) {              
          // for each person: find y-vector (list of targets $i$ has used)
          if (sender[n] != currentSender && sender[n] == x[m] && target[n] != currentTarget && attrvarCij[n] == attrCij ){
            y.push_back(target[n]);
          }
        } // n-loop close
        
        // clean up y (only unique values)
        std::sort( y.begin(), y.end() );
        y.erase( unique( y.begin(), y.end() ), y.end() );
        
        // interlock between x and y = wy vector
        wy.clear();
        std::sort( w.begin(), w.end() );
        std::sort( y.begin(), y.end() );
        std::set_intersection (w.begin(), w.end(), y.begin(), y.end(), std::back_inserter(wy) );
        // erase douplicates
        sort( wy.begin(), wy.end() );
        wy.erase( unique( wy.begin(), wy.end() ), wy.end() );
        
        // if these two actors (x-vector from a and y-vector from i) interact via a shared concept: calculate weightB (weightB = w_t(i, p))
        // positive/negative four cycle: only choose those events with different type
        if ( fourCycleType == "standard" ){
          if (wy.size() != 0) {
            for ( int o = begin; o < i; o++ ) {
              weightB = 0;
              if (sender[o] == x[m] && target[o] == currentTarget && attrvarBib[o] == attrBib && time[o] != currentTime) {
                weightB = std::abs(weightvar[o]) * exp( - ( currentTime - time[o] ) * xlog)  * xlog;
                tempTotalWeightB = tempTotalWeightB + weightB;
              }
            }// closes o-loop
          }
        }//closes if fourCycleType == standard
        if ( fourCycleType == "positive" ){
          if (wy.size() != 0) {
            for ( int o = begin; o < i; o++ ) {
              weightB = 0;
              if (sender[o] == x[m] && target[o] == currentTarget && typevar[o] == currentType && attrvarBib[o] == attrBib &&  time[o] != currentTime ) {
                weightB = std::abs(weightvar[o]) * exp( - ( currentTime - time[o] ) * xlog)  * xlog;
                tempTotalWeightB = tempTotalWeightB + weightB;
              }
            }// closes o-loop
          }
        }//closes if fourCycleType == positive
        if ( fourCycleType == "negative" ){
          if (wy.size() != 0) {
            for ( int o = begin; o < i; o++ ) {
              weightB = 0;
              if (sender[o] == x[m] && target[o] == currentTarget && typevar[o] != currentType && attrvarBib[o] == attrBib && time[o] != currentTime) {
                weightB = std::abs(weightvar[o]) * exp( - ( currentTime - time[o] ) * xlog)  * xlog;
                tempTotalWeightB = tempTotalWeightB + weightB;
              }
            }// closes o-loop
          }
        }//closes if fourCycleType == negative
        
        if ( fourCycleType == "standard"){
          // for each person: for each entry in wy: calculate weightA and weightC (p-loop open; q-loop open)
          for (int p = 0; p < wy.size(); p++) {
            
            tempTotalWeightA = 0;
            tempTotalWeightC = 0;
            
            for ( int q = begin; q < i; q++ ) {
              weightA = 0;
              weightC = 0;
              
              // calculate weightC (weightC = w_t(i, j))
              if (sender[q] == x[m] && target[q] == wy[p] && attrvarCij[q] == attrCij && time[q] != currentTime ) {
                weightC = std::abs(weightvar[q]) * exp( - ( currentTime - time[q] ) * xlog)  * xlog;
                tempTotalWeightC = tempTotalWeightC + weightC;
              }
              // calculate weightA (weightA = w_t(a, j))
              if (sender[q] == currentSender && target[q] == wy[p] && attrvarAaj[q] == attrAaj && time[q] != currentTime ) {
                weightA = std::abs(weightvar[q]) * exp( - ( currentTime - time[q] ) * xlog)  * xlog;
                tempTotalWeightA = tempTotalWeightA + weightA;
              }
            } //closes q-loop
            
            // close q-loop (backflash on each actor-concept-actor combination) & calculate weight
            // for each person (m) and concept in wy (p): 
              
              tempTotalWeightABC = tempTotalWeightA * tempTotalWeightB * tempTotalWeightC;
              
              // for each person: for each entry = sum up the multiplications
              totalWeightABC = tempTotalWeightABC + totalWeightABC;
          }//closes p-loop
        }//closes if fourCycleType == standard
        if ( fourCycleType == "positive" ){
          // for each person: for each entry in wy: calculate weightA and weightC (p-loop open; q-loop open)
          for (int p = 0; p < wy.size(); p++) {
            
            tempTotalWeightAPositive = 0;
            tempTotalWeightCPositive = 0;
            tempTotalWeightANegative = 0;
            tempTotalWeightCNegative = 0;
            
            for ( int q = begin; q < i; q++ ) {
              weightA = 0;
              weightC = 0;
              
              // if both weightA and weightC events are of same type, do this:
                // calculate weightC (weightC = w_t(i, j))
              if (sender[q] == x[m] && target[q] == wy[p] && typevar[q] == currentType && attrvarCij[q] == attrCij && time[q] != currentTime ) {
                weightC = std::abs(weightvar[q]) * exp( - ( currentTime - time[q] ) * xlog)  * xlog;
                tempTotalWeightCPositive = tempTotalWeightCPositive + weightC;
              }
              // calculate weightA (weightA = w_t(a, j))
              if (sender[q] == currentSender && target[q] == wy[p] && typevar[q] == currentType && attrvarAaj[q] == attrAaj && time[q] != currentTime ) {
                weightA = std::abs(weightvar[q]) * exp( - ( currentTime - time[q] ) * xlog)  * xlog;
                tempTotalWeightAPositive = tempTotalWeightAPositive + weightA;
              }
              // if both weightA and weightC events are of same type and negative, do this:
                // calculate weightC (weightC = w_t(i, j))
              if (sender[q] == x[m] && target[q] == wy[p] && typevar[q] != currentType && attrvarCij[q] == attrCij && time[q] != currentTime) {
                weightC = std::abs(weightvar[q]) * exp( - ( currentTime - time[q] ) * xlog)  * xlog;
                tempTotalWeightCNegative = tempTotalWeightCNegative + weightC;
              }
              // calculate weightA (weightA = w_t(a, j))
              if (sender[q] == currentSender && target[q] == wy[p] && typevar[q] != currentType && attrvarAaj[q] == attrAaj && time[q] != currentTime) {
                weightA = std::abs(weightvar[q]) * exp( - ( currentTime - time[q] ) * xlog) * xlog;
                tempTotalWeightANegative = tempTotalWeightANegative + weightA;
              }
            } //closes q-loop
            
            // close q-loop (backflash on each actor-concept-actor combination) & calculate weight
            // for each person (m) and concept in wy (p): 
              
              tempTotalWeightABCPositive = tempTotalWeightAPositive * tempTotalWeightB * tempTotalWeightCPositive;
              tempTotalWeightABCNegative = tempTotalWeightANegative * tempTotalWeightB * tempTotalWeightCNegative;
              
              // for each person: for each entry = sum up the multiplications
              totalWeightABCPositive = tempTotalWeightABCPositive + totalWeightABCPositive;
              totalWeightABCNegative = tempTotalWeightABCNegative + totalWeightABCNegative;
          }//closes p-loop
        }//closes if fourCycleType == positive
        if ( fourCycleType == "negative" ){
          // for each person: for each entry in wy: calculate weightA and weightC (p-loop open; q-loop open)
          for (int p = 0; p < wy.size(); p++) {
            
            tempTotalWeightAPositive = 0;
            tempTotalWeightCPositive = 0;
            tempTotalWeightANegative = 0;
            tempTotalWeightCNegative = 0;
            
            for ( int q = begin; q < i; q++ ) {
              weightA = 0;
              weightC = 0;
              
              // if both weightA and weightC events are of opposite type, do this:
                // calculate weightC (weightC = w_t(i, j))
              if (sender[q] == x[m] && target[q] == wy[p] && typevar[q] == currentType && attrvarCij[q] == attrCij && time[q] != currentTime ) {
                weightC = std::abs(weightvar[q]) * exp( - ( currentTime - time[q] ) * xlog)  * xlog;
                tempTotalWeightCPositive = tempTotalWeightCPositive + weightC;
              }
              // calculate weightA (weightA = w_t(a, j))
              if (sender[q] == currentSender && target[q] == wy[p] && typevar[q] != currentType && attrvarAaj[q] == attrAaj && time[q] != currentTime ) {
                weightA = std::abs(weightvar[q]) * exp( - ( currentTime - time[q] ) * xlog)  * xlog;
                tempTotalWeightAPositive = tempTotalWeightAPositive + weightA;
              }
              // if both weightA and weightC events are negative, do this:
                // calculate weightC (weightC = w_t(i, j))
              if (sender[q] == x[m] && target[q] == wy[p] && typevar[q] != currentType && attrvarCij[q] == attrCij && time[q] != currentTime) {
                weightC = std::abs(weightvar[q]) * exp( - ( currentTime - time[q] ) * xlog)  * xlog;
                tempTotalWeightCNegative = tempTotalWeightCNegative + weightC;
              }
              // calculate weightA (weightA = w_t(a, j))
              if (sender[q] == currentSender && target[q] == wy[p] && typevar[q] == currentType && attrvarAaj[q] == attrAaj && time[q] != currentTime ) {
                weightA = std::abs(weightvar[q]) * exp( - ( currentTime - time[q] ) * xlog) * xlog;
                tempTotalWeightANegative = tempTotalWeightANegative + weightA;
              }
            } //closes q-loop
            
            // close q-loop (backflash on each actor-concept-actor combination) & calculate weight
            // for each person (m) and concept in wy (p): 
              
              tempTotalWeightABCPositive = tempTotalWeightAPositive * tempTotalWeightB * tempTotalWeightCPositive;
              tempTotalWeightABCNegative = tempTotalWeightANegative * tempTotalWeightB * tempTotalWeightCNegative;
              
              // for each person: for each entry = sum up the multiplications
              totalWeightABCPositive = tempTotalWeightABCPositive + totalWeightABCPositive;
              totalWeightABCNegative = tempTotalWeightABCNegative + totalWeightABCNegative;
          }//closes p-loop
        }//closes if fourCycleType == negative
      } // m-loop close
      
      if ( fourCycleType == "standard"){
        totalWeightABC = std::pow(totalWeightABC, 1/3.); //whyever, there has to be a . behind the 1/3
      }else{
        totalWeightABC = std::pow(totalWeightABCPositive, 1/3.) + std::pow(totalWeightABCNegative, 1/3.); //wieso auch immer - aber da muss ein Punkt hinter die 3
      }
      result = totalWeightABC; 
     
  // attrvarNow => in R-loop now 
  //}else{ //closes "if ( attrvarNow[i] == attrNow ) {}"
  //  result = 0.0;
  //}
  return result;
  }


//####################################################################
// [[Rcpp::export]]
double similarityTotalAverageCpp(
  std::vector<std::string> sender,
  std::string currentSender,
  std::vector<std::string> target,
  std::string currentTarget,
  NumericVector time,
  double currentTime,
  std::vector<std::string> eventAttributeVar,
  std::string eventAttribute,
  std::vector<std::string> eventTypeVar,
  std::string currentType,
  std::string totalAverageSim,
  std::string matchNomatchSim,
  std::string senderTargetSim, 
  std::vector<std::string> v, // sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
  std::vector<std::string> w, // sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
  int i, 
  int begin) {
  
  double result = 0.0;
  std::vector<std::string> x;
  std::vector<std::string> xw;
  std::vector<double> a_positive;
  std::vector<double> a_negative;
  std::vector<double> i_negative;
  std::vector<double> i_positive;	
  std::vector<std::string> xwneg;
  double totalNumber;
  double numberNoMatch;
  double numberMatch;
    
    i_negative.clear();
    a_negative.clear();
    i_positive.clear();
    a_positive.clear();
    xwneg.clear();
    totalNumber = 0;
       
    // for each entry in v
    for (int k = 0; k < v.size(); k++){
      
      x.clear();
      xw.clear();
      
      if ( senderTargetSim == "sender" ){
        for (int l = begin; l < i; l++) {
          // if the event has concept k in v
          if (sender[l] == v[k] && time[l] != currentTime && target[l] != currentTarget && eventAttributeVar[l] == eventAttribute)  {      
            x.push_back(target[l]);
          }
        }//closes l-loop
      }else{
        for (int l = begin; l < i; l++) {
          // if the event has concept k in v
          if (target[l] == v[k] && time[l] != currentTime && sender[l] != currentSender && eventAttributeVar[l] == eventAttribute)  {      
            x.push_back(sender[l]);
          }
        }//closes l-loop
      }
      
      // clean up x-vector
      std::sort( x.begin(), x.end() );
      x.erase( unique( x.begin(), x.end() ), x.end() );
      
      // for each entry in v => get intersection between x and w (w = actors who said p/what else has actor a said?) = filter actors who said p and v[k]
      xw.clear();
      std::sort( w.begin(), w.end() );
      std::sort( x.begin(), x.end() );
      std::set_intersection (x.begin(), x.end(), w.begin(), w.end(), std::back_inserter(xw) );
      // erase douplicates
      sort( xw.begin(), xw.end() );
      xw.erase( unique( xw.begin(), xw.end() ), xw.end() );
      
      if ( matchNomatchSim == "match"){
        
        // if there acctually is an intersection in xw
        if (xw.size() != 0 ) {
          
          // for each entry in xw:
            for (int m = 0; m < xw.size(); m++) {
              
              i_negative.clear();
              i_positive.clear();
              a_negative.clear();
              a_positive.clear();
              
              // loop back over all events until i
              for (int n = begin; n < i; n++){
                if ( senderTargetSim == "sender" ){
                  if (sender[n] == v[k] && target[n] == xw[m] && eventTypeVar[n] == currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    i_positive.push_back(time[n]);
                  }
                  if (sender[n] == v[k] && target[n] == xw[m] && eventTypeVar[n] != currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    i_negative.push_back(time[n]);
                  }
                  if (sender[n] == currentSender && target[n] == xw[m] && eventTypeVar[n] == currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    a_positive.push_back(time[n]);
                  }
                  if (sender[n] == currentSender && target[n] == xw[m] && eventTypeVar[n] != currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    a_negative.push_back(time[n]);
                  }
                }else{
                  if (target[n] == v[k] && sender[n] == xw[m] && eventTypeVar[n] == currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    i_positive.push_back(time[n]);
                  }
                  if (target[n] == v[k] && sender[n] == xw[m] && eventTypeVar[n] != currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    i_negative.push_back(time[n]);
                  }
                  if (target[n] == currentTarget && sender[n] == xw[m] && eventTypeVar[n] == currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    a_positive.push_back(time[n]);
                  }
                  if (target[n] == currentTarget && sender[n] == xw[m] && eventTypeVar[n] != currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    a_negative.push_back(time[n]);
                  }
                }//closes senderTargetSim == target
              }//closes n-loop		
            }//closes m-loop
          
        }//closes if xw.size != 0
        
        numberNoMatch = 0;
        numberMatch = 0;
        
        //how large are the respextive vectors with the matches in them?
        if (a_positive.size() >= i_positive.size() && i_positive.size() != 0 ) {
          numberMatch = i_positive.size();
        }
        if (a_negative.size() >= i_negative.size() && i_negative.size() != 0) {
          numberNoMatch = i_negative.size();
        }
        if (i_positive.size() > a_positive.size() && a_positive.size() != 0) {
          numberMatch = a_positive.size();
        }
        if (i_negative.size() > a_negative.size() && a_negative.size() != 0) {
          numberNoMatch = a_negative.size();
        }
        //// how many actors used concept v[k] in same manner as a? // how many concepts are used by both in same manner?
        totalNumber = totalNumber + numberNoMatch + numberMatch;
      }else{ // if matchNomatchSim = "nomatch"
             // how many actors used concept v[k]? // how many concepts are used by both?
             totalNumber = totalNumber + xw.size();
      }
      
    }//closes k-loop
    
      if (totalAverageSim == "total") {
        result = totalNumber;
      }
      if (totalAverageSim == "average") {
        result = totalNumber/v.size(); //TODO: correct to divide by v.size?
      }
    
  return result;
  }


//####################################################################
// [[Rcpp::export]]
double similaritySimpleCpp(
  std::vector<std::string> sender,
  std::string currentSender,
  std::vector<std::string> target,
  std::string currentTarget,
  NumericVector time,
  double currentTime,
  double xlog,
  std::vector<std::string> eventAttributeVar,
  std::string eventAttribute,
  std::vector<std::string> eventTypeVar,
  std::string currentType,
  std::string matchNomatchSim,
  std::string senderTargetSim, 
  std::vector<std::string> v, // sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
  std::vector<std::string> w, // sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
  int i, 
  int begin) {
  
  double result = 0;
  std::vector<std::string> x;
  std::vector<std::string> xw;
  double a_positive;
  double a_negative;
  double i_negative;
  double i_positive;	
  std::vector<std::string> xwneg;
  double totalNumber = 0;
  double timePLast = 0;
  int counter = 0;
  double totalSim = 0;
  double weightSim;
  
  // for each entry in v
  for (int k = 0; k < v.size(); k++){
    
    x.clear();
    xw.clear();
    totalNumber = 0;
    
    if ( senderTargetSim == "sender" ){
      for (int l = begin; l < i; l++) {
        // if the event has concept k in v
        if (sender[l] == v[k] && time[l] != currentTime && target[l] != currentTarget && eventAttributeVar[l] == eventAttribute)  {      
          x.push_back(target[l]);
        }
      }//closes l-loop
    }else{
      for (int l = begin; l < i; l++) {
        // if the event has concept k in v
        if (target[l] == v[k] && time[l] != currentTime && sender[l] != currentSender && eventAttributeVar[l] == eventAttribute)  {      
          x.push_back(sender[l]);
        }
      }//closes l-loop
    }
    
    // clean up x-vector
    std::sort( x.begin(), x.end() );
    x.erase( unique( x.begin(), x.end() ), x.end() );
    
    // for each entry in v => get intersection between x and w (w = actors who said p/what else has actor a said?) = filter actors who said p and v[k]
    xw.clear();
    std::sort( w.begin(), w.end() );
    std::sort( x.begin(), x.end() );
    std::set_intersection (x.begin(), x.end(), w.begin(), w.end(), std::back_inserter(xw) );
    // erase douplicates
    sort( xw.begin(), xw.end() );
    xw.erase( unique( xw.begin(), xw.end() ), xw.end() );
    
    
    // Match: check for each overlaping actor/target in xw, whether the type matches
    if ( matchNomatchSim == "match"){
      
      // if there acctually is an intersection in xw
      if (xw.size() != 0 ) {
        
        // for each entry in xw:
          for (int m = 0; m < xw.size(); m++) {
            //
              i_negative = 0;
              i_positive = 0;
              a_negative = 0;
              a_positive = 0;
              
              //check if a and i used the same type. If they used a certan type = give them a 1 - then compare if both a  and i have 1 in the same type-cateogry
              for (int n = begin; n < i; n++){
                if ( senderTargetSim == "sender" ){
                  if (sender[n] == v[k] && target[n] == xw[m] && eventTypeVar[n] == currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    i_positive = 1;
                  }
                  if (sender[n] == v[k] && target[n] == xw[m] && eventTypeVar[n] != currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    i_negative = 1; 
                  }
                  if (sender[n] == currentSender && target[n] == xw[m] && eventTypeVar[n] == currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    a_positive = 1;
                  }
                  if (sender[n] == currentSender && target[n] == xw[m] && eventTypeVar[n] != currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    a_negative = 1;
                  }
                }else{
                  if (target[n] == v[k] && sender[n] == xw[m] && eventTypeVar[n] == currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    i_positive = 1;
                  }
                  if (target[n] == v[k] && sender[n] == xw[m] && eventTypeVar[n] != currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    i_negative = 1;
                  }
                  if (target[n] == currentTarget && sender[n] == xw[m] && eventTypeVar[n] == currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    a_positive = 1; 
                  }
                  if (target[n] == currentTarget && sender[n] == xw[m] && eventTypeVar[n] != currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    a_negative = 1;
                  }
                }//closes senderTargetSim == target
              }//closes n-loop
              // compare a_positive with i_positive and a_negative with i_negative, if they match, add +1 to totalNumber
              // technically: totalNumber can be twice xw.size() => if all actors/targets use both positive and negative targets
              if(a_positive == 1 && i_positive == 1){
                totalNumber = totalNumber + 1;
              }
              if(a_negative == 1 && i_negative == 1){
                totalNumber = totalNumber + 1;
              }
          }//closes m-loop
      }//closes if xw.size != 0
      
    }else{ // if matchNomatchSim = "nomatch"
    
    // how many actors used concept v[k]? // how many concepts are used by both?
    totalNumber = totalNumber + xw.size();
    }//closes matchNomatchSim == "nomatch"
    
    // how many actors used concept v[k]? // how many concepts are used by both?
    if ( xw.size() != 0 ){
      // find time for the time-discount in the simple-similarity equation
      if ( senderTargetSim == "sender"){
        // find time, when actor k used concept $p$ last
        for (int q = i-1; q >= begin; q--){
          if ( matchNomatchSim == "nomatch"){
            if (sender[q] == v[k] && target[q] == currentTarget && time[q] != currentTime &&
                eventAttributeVar[q] == eventAttribute){
              timePLast = time[q];
              break;
            }
          }else if (matchNomatchSim == "match"){
            if (sender[q] == v[k] && target[q] == currentTarget && time[q] != currentTime &&
                eventAttributeVar[q] == eventAttribute && eventTypeVar[q] == currentType ){
              timePLast = time[q];
              break;
            }
          }
        }//closes q-loop
      }else if ( senderTargetSim == "target"){
        // find time, when target k was last used by $a$
          for (int q = i-1; q >= begin; q--){
            if (matchNomatchSim == "nomatch"){
              if (target[q] == v[k] && sender[q] == currentSender && time[q] != currentTime &&
                  eventAttributeVar[q] == eventAttribute){
                timePLast = time[q];
                break;
              }
            } else if (matchNomatchSim == "match"){
              if (target[q] == v[k] && sender[q] == currentSender && time[q] != currentTime &&
                  eventAttributeVar[q] == eventAttribute && eventTypeVar[q] == currentType ){
                timePLast = time[q];
                break;
              }
            }
          }//closes q-loop
      }//closes senderTargetSim == "target"
      
      // calculate weight of each actor/target
      weightSim = totalNumber * exp(-(currentTime - timePLast)*xlog) * xlog;
      totalSim = totalSim + weightSim;
      if (weightSim != 0) {
        counter++;   
      }
    }//closes if xw.size != 0
  }//closes k-loop
  
  if (counter == 0) {
    result = 0.0;
  }else{
    result = (totalSim/counter);
  }
  return result;
  }


//####################################################################
// [[Rcpp::export]]
double similarityComplexCpp(
  std::vector<std::string> sender,
  std::string currentSender,
  std::vector<std::string> target,
  std::string currentTarget,
  NumericVector time,
  double currentTime,
  double xlog,
  double halflifeTimeDifference,
  std::vector<std::string> eventAttributeVar,
  std::string eventAttribute,
  std::vector<std::string> eventTypeVar,
  std::string currentType,
  std::string matchNomatchSim,
  std::string senderTargetSim,
  std::vector<std::string> v, // sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
  std::vector<std::string> w, // sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
  int i, 
  int begin) {
  
  double result = 0.0;
  std::vector<std::string> x;
  std::vector<std::string> xw;
  std::vector<std::string> xwneg;
  std::vector<double> a_positive;
  std::vector<double> a_negative;
  std::vector<double> i_negative;
  std::vector<double> i_positive;	
  std::vector<double> i_sendertarget;
  std::vector<double> a_sendertarget;
  double timePLast = 0.0;
  int counter = 0;
  double totalSim = 0.0;
  double weightSim;
  double sumCouplePositive;
  double sumCoupleNegative;
  double couple;
  double sumCoupleNoMatch;
  double sumConcept = 0.0;
      
    // for each entry in v
    for (int k = 0; k < v.size(); k++){
      
      x.clear();
      xw.clear();
      
      if ( senderTargetSim == "sender" ){
        for (int l = 0; l < i; l++) {
          // if the event has concept k in v
          if (sender[l] == v[k] && time[l] != currentTime && target[l] != currentTarget && eventAttributeVar[l] == eventAttribute)  {      
            x.push_back(target[l]);
          }
        }//closes l-loop
      }else{
        for (int l = 0; l < i; l++) {
          // if the event has concept k in v
          if (target[l] == v[k] && time[l] != currentTime && sender[l] != currentSender && eventAttributeVar[l] == eventAttribute)  {      
            x.push_back(sender[l]);
          }
        }//closes l-loop
      }
      
      // clean up x-vector
      std::sort( x.begin(), x.end() );
      x.erase( unique( x.begin(), x.end() ), x.end() );
      
      // for each entry in v => get intersection between x and w (w = actors who said p/what else has actor a said?) = filter actors who said p and v[k]
      xw.clear();
      std::sort( w.begin(), w.end() );
      std::sort( x.begin(), x.end() );
      std::set_intersection (x.begin(), x.end(), w.begin(), w.end(), std::back_inserter(xw) );
      // erase douplicates
      sort( xw.begin(), xw.end() );
      xw.erase( unique( xw.begin(), xw.end() ), xw.end() );
      
      if ( matchNomatchSim == "match"){
        
        // if there acctually is an intersection in xw
        if (xw.size() != 0 ) {
          
          // for each entry in xw:
            for (int m = 0; m < xw.size(); m++) {
              
              i_negative.clear();
              i_positive.clear();
              a_negative.clear();
              a_positive.clear();
              i_sendertarget.clear();
              a_sendertarget.clear();
              
              // loop back over all events until i
              for (int n = 0; n < i; n++){
                if ( senderTargetSim == "sender" ){
                  if (sender[n] == v[k] && target[n] == xw[m] && eventTypeVar[n] == currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    i_positive.push_back(time[n]);
                  }
                  if (sender[n] == v[k] && target[n] == xw[m] && eventTypeVar[n] != currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    i_negative.push_back(time[n]);
                  }
                  if (sender[n] == currentSender && target[n] == xw[m] && eventTypeVar[n] == currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    a_positive.push_back(time[n]);
                  }
                  if (sender[n] == currentSender && target[n] == xw[m] && eventTypeVar[n] != currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    a_negative.push_back(time[n]);
                  }
                }else{
                  if (target[n] == v[k] && sender[n] == xw[m] && eventTypeVar[n] == currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    i_positive.push_back(time[n]);
                  }
                  if (target[n] == v[k] && sender[n] == xw[m] && eventTypeVar[n] != currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    i_negative.push_back(time[n]);
                  }
                  if (target[n] == currentTarget && sender[n] == xw[m] && eventTypeVar[n] == currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    a_positive.push_back(time[n]);
                  }
                  if (target[n] == currentTarget && sender[n] == xw[m] && eventTypeVar[n] != currentType && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                    a_negative.push_back(time[n]);
                  }
                }//closes senderTargetSim == target
              }//closes n-loop	
            }//closes m-loop
        }//closes if xw.size != 0
        
        sumCouplePositive = 0;
        sumCoupleNegative = 0;
        couple = 0;
        
        //how large are the respextive vectors with the matches in them?
        if (a_positive.size() >= i_positive.size() && i_positive.size() != 0 ) {
          for (int p = 0; p < i_positive.size(); p++){
            couple = 1 * exp(-(std::abs(i_positive[p]-a_positive[p]))*(log(2.0)/halflifeTimeDifference));
            sumCouplePositive = sumCouplePositive + couple;
            couple = 0;
          }
        }
        if (a_negative.size() >= i_negative.size() && i_negative.size() != 0) {
          for (int p = 0; p < i_negative.size(); p++){
            couple = 1 * exp(-(std::abs(i_negative[p]-a_negative[p]))*(log(2.0)/halflifeTimeDifference));
            sumCoupleNegative = sumCoupleNegative + couple;
            couple = 0;
          }
        }
        if (i_positive.size() > a_positive.size() && a_positive.size() != 0) {
          for (int p = 0; p < a_positive.size(); p++){
            couple = 1 * exp(-(std::abs(i_positive[p]-a_positive[p]))*(log(2.0)/halflifeTimeDifference));
            sumCouplePositive = sumCouplePositive + couple;
            couple = 0;
          }
        }
        if (i_negative.size() > a_negative.size() && a_negative.size() != 0) {
          for (int p = 0; p < a_negative.size(); p++){
            couple = 1 * exp(-(std::abs(i_negative[p]-a_negative[p]))*(log(2.0)/halflifeTimeDifference));
            sumCoupleNegative = sumCoupleNegative + couple;
            couple = 0;
          }
        }
        //// how many actors used concept v[k] in same manner as a? // how many concepts are used by both in same manner?
        sumConcept = sumConcept + sumCoupleNegative + sumCouplePositive;
      }else{ // if matchNomatchSim = "nomatch"
             // if there acctually is an intersection in xw
             if (xw.size() != 0 ) {
               
               // for each entry in xw:
                 for (int m = 0; m < xw.size(); m++) {
                   i_sendertarget.clear();
                   a_sendertarget.clear();
                   
                   // loop back over all events until i
                   for (int n = 0; n < i; n++){
                     if ( senderTargetSim == "sender" ){
                       if (sender[n] == v[k] && target[n] == xw[m] && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                         i_sendertarget.push_back(time[n]);
                       }
                       if (sender[n] == currentSender && target[n] == xw[m]  && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                         a_sendertarget.push_back(time[n]);
                       }
                     }else{
                       if (target[n] == v[k] && sender[n] == xw[m] && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                         i_sendertarget.push_back(time[n]);
                       }
                       if (target[n] == currentTarget && sender[n] == xw[m] && time[n] != currentTime && eventAttributeVar[n] == eventAttribute){
                         a_sendertarget.push_back(time[n]);
                       }
                     }//closes senderTargetSim == target
                   }//closes n-loop
                 }//closes m-loop
             }//closes xw.size() != 0
             
             sumCoupleNoMatch = 0;
             
             //how large are the respextive vectors with the matches in them?
             if (a_sendertarget.size() >= i_sendertarget.size() && i_sendertarget.size() != 0) {
               for (int p = 0; p < i_negative.size(); p++){
                 couple = 1 * exp(-(std::abs(i_sendertarget[p]-a_sendertarget[p]))*(log(2.0)/halflifeTimeDifference));
                 sumCoupleNoMatch = sumCoupleNoMatch + couple;
                 couple = 0;
               }
             }
             if (i_sendertarget.size() > a_sendertarget.size() && a_sendertarget.size() != 0) {
               for (int p = 0; p < a_positive.size(); p++){
                 couple = 1 * exp(-(std::abs(i_sendertarget[p]-a_sendertarget[p]))*(log(2.0)/halflifeTimeDifference));
                 sumCoupleNoMatch = sumCoupleNoMatch + couple;
                 couple = 0;
               }
             }
             //add them together
             sumConcept = sumConcept + sumCoupleNoMatch;	
      }//closes matchNomatchSim == "nomatch"
      
      // how many actors used concept v[k]? // how many concepts are used by both?
      if ( xw.size() != 0 ){
        // find time for the time-discount in the simple-similarity equation
        if ( senderTargetSim == "sender"){
          // find time, when actor k used concept $p$ last
          for (int q = i-1; q >= 0; q--){
            if ( matchNomatchSim == "nomatch"){
              if (sender[q] == v[k] && target[q] == currentTarget && time[q] != currentTime &&
                    eventAttributeVar[q] == eventAttribute){
                timePLast = time[q];
                break;
              }
            }else if (matchNomatchSim == "match"){
              if (sender[q] == v[k] && target[q] == currentTarget && time[q] != currentTime &&
                    eventAttributeVar[q] == eventAttribute && eventTypeVar[q] == currentType){
                timePLast = time[q];
                break;
              }
            }
          }//closes q-loop
        }else if ( senderTargetSim == "target"){
          // find time, when target k was last used by $a$
            for (int q = i-1; q >= 0; q--){
              if (matchNomatchSim == "nomatch"){
                if (target[q] == v[k] && sender[q] == currentSender && time[q] != currentTime &&
                      eventAttributeVar[q] == eventAttribute){
                  timePLast = time[q];
                  break;
                }
              } else if (matchNomatchSim == "match"){
                if (target[q] == v[k] && sender[q] == currentSender && time[q] != currentTime &&
                      eventAttributeVar[q] == eventAttribute && eventTypeVar[q] == currentType){
                  timePLast = time[q];
                  break;
                }
              }
            }//closes q-loop
        }//closes senderTargetSim == "target"
        
        // calculate weight of each actor/target
        weightSim = sumConcept * exp(-(currentTime-timePLast)*xlog) * xlog;
        totalSim = totalSim + weightSim;
        if (weightSim != 0) {
          counter++;   
        }
      }//closes if xw.size != 0
    }//closes k-loop
    
    if (counter == 0) {
      result = 0.0;
    }else{
      result = (totalSim/counter);
    }
  return result;
  }


//####################################################################
// [[Rcpp::export]]
NumericVector reciprocityCpp(
  NumericVector time,
  NumericVector weightvar,
  std::vector<std::string> sender,
  std::vector<std::string> target,
  std::vector<std::string> typevar,
  std::string type1, 
  std::string type2, 
  std::vector<std::string> attrvar,
  std::string attr1,
  std::string attr2, 
  double xlog, 
  std::string reciprocitytype ) {
  
  NumericVector result(sender.size());
  
  // for-loop i: for each event, do:
  for ( int i = 0; i < sender.size(); i++){
      
      // reset all the variables
      double reciprocity = 0;
      double totalreciprocity = 0;
      double weight = 0;
      double resulttemp = 0;
      
      // for-loop w: go back over all events and filter them
      for ( int w = 0; w < i; w++ ){
        
        // sender-target reciprocity only
        if ( reciprocitytype == "s-t-only" ){
          if ( target[w] == sender[i] && sender[w] == target[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typematch 
        if ( reciprocitytype == "s-t-typematch" ){
          if ( target[w] == sender[i] && sender[w] == target[i] && typevar[w] == typevar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typemix
        if ( reciprocitytype == "s-t-typemix" ){
          if ( target[w] == sender[i] && sender[w] == target[i] && typevar[i] == type1 && typevar[w] == type2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typefilter
        if ( reciprocitytype == "s-t-typefilter" ){
          if ( target[w] == sender[i] && sender[w] == target[i] && typevar[w] == type1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-attributematch
        if ( reciprocitytype == "s-t-attributematch" ){
          if ( target[w] == target[i] && sender[w] == sender[i] && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-attributemix
        if ( reciprocitytype == "s-t-attributemix" ){
          if ( target[w] == sender[i] && sender[w] == target[i] && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-attributefilter
        if ( reciprocitytype == "s-t-attributefilter" ){
          if ( target[w] == sender[i] && sender[w] == target[i] && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typematch-attributematch
        if ( reciprocitytype == "s-t-typematch-attributematch" ){
          if ( target[w] == sender[i] && sender[w] == target[i] && typevar[w] == typevar[i] && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typematch-attributemix
        if ( reciprocitytype == "s-t-typematch-attributemix" ){
          if ( target[w] == sender[i] && sender[w] == target[i] && typevar[w] == typevar[i] && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typematch-attributefilter
        if ( reciprocitytype == "s-t-typematch-attributefilter" ){
          if ( target[w] == sender[i] && sender[w] == target[i] && typevar[w] == typevar[i] && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typemix-attributematch
        if ( reciprocitytype == "s-t-typemix-attributematch" ){
          if ( target[w] == sender[i] && sender[w] == target[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typemix-attributemix
        if ( reciprocitytype == "s-t-typemix-attributemix" ){
          if ( target[w] == sender[i] && sender[w] == target[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typemix-attributefilter
        if ( reciprocitytype == "s-t-typemix-attributefilter" ){
          if ( target[w] == sender[i] && sender[w] == target[i] && typevar[i] == type1 && typevar[w] == type2 && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typefilter-attributematch
        if ( reciprocitytype == "s-t-typefilter-attributematch" ){
          if ( target[w] == sender[i] && sender[w] == target[i] && typevar[w] == type1 && attrvar[w] == attrvar[i] ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typefilter-attributemix
        if ( reciprocitytype == "s-t-typefilter-attributemix" ){
          if ( target[w] == sender[i] && sender[w] == target[i] && typevar[w] == type1 && attrvar[i] == attr1 && attrvar[w] == attr2 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        // sender-target-typefilter-attributefilter
        if ( reciprocitytype == "s-t-typefilter-attributefilter" ){
          if ( target[w] == sender[i] && sender[w] == target[i] && typevar[w] == type1 && attrvar[w] == attr1 ){
            weight = std::abs(weightvar[w]);
          } else {
            weight = 0;
          }
        }
        
        //
          // calculate the weight
        if ( time[i] == time[w] ){
          totalreciprocity = totalreciprocity + 0;
        } else {
          reciprocity = weight * exp( - ( time[i] - time[w] ) * xlog) * xlog ;
          totalreciprocity = totalreciprocity + reciprocity;
        }
        
        resulttemp = totalreciprocity;
      } // closes w-loop
      
      // ascribe calculated weights to the result-variable  	
      result[i] = resulttemp;
      
    } // closes i-loop
  
  // return variable and hand it back to R
  return Rcpp::wrap(result);
}


//####################################################################
// [[Rcpp::export]]
double triadCpp(
	std::vector<std::string> v, // intersection x and y
	std::vector<std::string> sender, 
	std::vector<std::string> target,
	NumericVector time, 
	NumericVector weightvar, 
	std::vector<std::string> typevar, 
	std::string typeA, 
	std::string typeB,
	std::vector<std::string> attributevarAI,
	std::string attrAI,
	std::vector<std::string> attributevarBI,
	std::string attrBI,
    double xlog, 
	int i, 
	std::string currentSender, 
	std::string currentTarget, 
	double currentTime) {
		
	double result = 0.0;
	double weighta;
	double weightb;
	double totalweighta;
	double totalweightb;
	double weightab = 0.0;
	double totalweight = 0.0;
	
	// for each entry in v
	for (int j = 0; j < v.size(); j++) {
		totalweighta = 0.0;
		totalweightb = 0.0;
		for ( int z = 0; z < i-1; z++ ) {   
			weighta = 0.0;
			weightb = 0.0;
			//caluculate weighta
			if ( ( (sender[z] == currentSender && target[z] == v[j]) || 
				(target[z] == currentSender && sender[z] == v[j]) ) && typevar[z] == typeA &&
			attributevarAI[z] == attrAI && time[z] != currentTime ){
				weighta = std::abs(weightvar[z]) * exp( - ( currentTime - time[z] ) * xlog) * xlog;
				totalweighta = totalweighta + weighta;   
			}
			//calculate weightb
			if ( ( (sender[z] == currentTarget && target[z] == v[j]) || 
				(target[z] == currentTarget && sender[z] == v[j]) ) && typevar[z] == typeB &&
			attributevarBI[z] == attrBI && time[z] != currentTime ){
				weightb = std::abs(weightvar[z]) * exp( - ( currentTime- time[z] ) * xlog) * xlog;
				totalweightb = totalweightb + weightb;   
			}
		} //closes z-loop

		// multiply totalweighta times totalweightb
		weightab = totalweighta * totalweightb;
		totalweight = totalweight + weightab;   
	} //closes j-loop  
	//take the squared rood of totalweight
	result = sqrt(totalweight);
	return result;
  }


//####################################################################
// [[Rcpp::export]]
double weightTimesSummationCpp(
  NumericVector pastSenderTimes,
  double xlog, 
  double currentTime, 
  NumericVector weightvar) {

  double totalWeight = 0.0;
  double weight = 0.0;
  double result = 0.0;

  // for each bill that the current sender has cosponsored in past
  for (int j = 0; j < pastSenderTimes.size(); j++){
    weight = weightvar[j] * exp( - ( currentTime - pastSenderTimes[j] ) * xlog)  * xlog;
    totalWeight = totalWeight + weight ;
  }// close j-loop

  result = totalWeight;
  return result;
}


//####################################################################
// [[Rcpp::export]]
DataFrame createNullEvents(
    std::vector<std::string> eventID,
    std::vector<std::string> sender,
    std::vector<std::string> target,
    std::vector<std::string> eventAttribute,
    std::vector<double> time,
    std::vector<double> start,
    std::vector<double> end, 
    std::vector<double> allEventTimes) {
  
  DataFrame result;
  std::vector<std::string> eventIDNew;
  std::vector<std::string> senderNew;
  std::vector<std::string> targetNew;
  std::vector<std::string> eventAttributeNew;
  NumericVector startNew;
  NumericVector endNew;
  NumericVector eventTime;
  NumericVector eventDummy;
  
  //for each event in the sequence
  for (int i = 0; i < sender.size(); i++){
    
    //for each event in allEventTimes
    for(int w = 0; w < allEventTimes.size(); w++){
      
      //for each eventtime between start[i] and end[i]
      if(allEventTimes[w] >= start[i] && allEventTimes[w] <= end[i]){
        
        eventIDNew.push_back(eventID[i]);
        senderNew.push_back(sender[i]);
        targetNew.push_back(target[i]);
        eventAttributeNew.push_back(eventAttribute[i]);
        startNew.push_back(start[i]);
        endNew.push_back(end[i]);
        eventTime.push_back(allEventTimes[w]);
        //is it a null-event or a true-event?
        if(time[i] == allEventTimes[w]){
          eventDummy.push_back(1);
        }else{
          eventDummy.push_back(0);
        }
        
      }
    }//closes w-loop
  }//closes i-loop
  
  //combine all vectors into one
  result = Rcpp::DataFrame::create(Rcpp::Named("eventID") = eventIDNew, 
                                   Rcpp::Named("sender") = senderNew, 
                                   Rcpp::Named("target") = targetNew,
                                   Rcpp::Named("eventTime") = eventTime,
                                   Rcpp::Named("eventDummy") = eventDummy, 
                                   Rcpp::Named("eventAtRiskFrom") = startNew, 
                                   Rcpp::Named("eventAtRiskUntil") = endNew, 
                                   Rcpp::Named("eventAttribute") = eventAttributeNew);
  return Rcpp::wrap(result);
}