#include <string.h>
#include <Rcpp.h>
using namespace Rcpp;


//####################################################################
// [[Rcpp::export]]
NumericVector absoluteDiffAverageWeightEventAttributeCpp(
    std::vector<std::string> sender,
    std::vector<std::string> target,
    NumericVector time,
    NumericVector weightvar,
    std::vector<std::string> eventattributevar,
    std::string eventattribute,
    double xlog) {
    
    NumericVector result(time.size());
    int count = 0;
    
    // for-loop: for each event, do:
      for ( int i = 0; i < time.size(); i++){
        
        double totaldiff = 0;
        count = 0;
        
        for ( int w = 0; w < i; w++ ) {
          if (eventattributevar[w] == eventattribute && sender[i] != sender[w] && target[i] == target[w]) {
            if ( time[i] != time[w] ) {
              count = count + 1;
              totaldiff = totaldiff + std::abs(weightvar[i] - weightvar[w]);
            } 
          }
        } // closes w-loop
        if (count == 0) {
          result[i] = 0;
        }else{
          result[i] = totaldiff/count;
        }
      } // closes i-loop
    return Rcpp::wrap(result);
}













