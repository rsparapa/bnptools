hbart_vartivity = function(
                              ip, ##number of predictor variables
                              im,  ##number of trees in mean model
                              imh, ##number of trees in variance model
                              ind, ##number of draws saved from the posterior
                              itc, ##number of parallel compute threads
                              ifit ##saved fitted model object returned from cpsambrt
                              )
{
    res=.Call("cpsambrt_vartivity",
              ip, ##number of predictor variables
              im,  ##number of trees in mean model
              imh, ##number of trees in variance model
              ind, ##number of draws saved from the posterior
              itc, ##number of parallel compute threads
              ifit, ##saved fitted model object returned from cpsambrt
              PACKAGE='hbart')
    return(res)
}
