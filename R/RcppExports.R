# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

fourCycleCpp <- function(sender, currentSender, target, currentTarget, typevar, currentType, time, currentTime, weightvar, xlog, attrvarAaj, attrAaj, attrvarBib, attrBib, attrvarCij, attrCij, fourCycleType, w, x, i, begin) {
    .Call(`_rem_fourCycleCpp`, sender, currentSender, target, currentTarget, typevar, currentType, time, currentTime, weightvar, xlog, attrvarAaj, attrAaj, attrvarBib, attrBib, attrvarCij, attrCij, fourCycleType, w, x, i, begin)
}

similarityTotalAverageCpp <- function(sender, currentSender, target, currentTarget, time, currentTime, eventAttributeVar, eventAttribute, eventTypeVar, currentType, totalAverageSim, matchNomatchSim, senderTargetSim, v, w, i, begin) {
    .Call(`_rem_similarityTotalAverageCpp`, sender, currentSender, target, currentTarget, time, currentTime, eventAttributeVar, eventAttribute, eventTypeVar, currentType, totalAverageSim, matchNomatchSim, senderTargetSim, v, w, i, begin)
}

similaritySimpleCpp <- function(sender, currentSender, target, currentTarget, time, currentTime, xlog, eventAttributeVar, eventAttribute, eventTypeVar, currentType, matchNomatchSim, senderTargetSim, v, w, i, begin) {
    .Call(`_rem_similaritySimpleCpp`, sender, currentSender, target, currentTarget, time, currentTime, xlog, eventAttributeVar, eventAttribute, eventTypeVar, currentType, matchNomatchSim, senderTargetSim, v, w, i, begin)
}

similarityComplexCpp <- function(sender, currentSender, target, currentTarget, time, currentTime, xlog, halflifeTimeDifference, eventAttributeVar, eventAttribute, eventTypeVar, currentType, matchNomatchSim, senderTargetSim, v, w, i, begin) {
    .Call(`_rem_similarityComplexCpp`, sender, currentSender, target, currentTarget, time, currentTime, xlog, halflifeTimeDifference, eventAttributeVar, eventAttribute, eventTypeVar, currentType, matchNomatchSim, senderTargetSim, v, w, i, begin)
}

triadCpp <- function(v, sender, target, time, weightvar, typevar, typeA, typeB, attributevarAI, attrAI, attributevarBI, attrBI, xlog, i, currentSender, currentTarget, currentTime) {
    .Call(`_rem_triadCpp`, v, sender, target, time, weightvar, typevar, typeA, typeB, attributevarAI, attrAI, attributevarBI, attrBI, xlog, i, currentSender, currentTarget, currentTime)
}

weightTimesSummationCpp <- function(pastSenderTimes, xlog, currentTime, weightvar) {
    .Call(`_rem_weightTimesSummationCpp`, pastSenderTimes, xlog, currentTime, weightvar)
}

createNullEvents <- function(eventID, sender, target, eventAttribute, time, start, end, allEventTimes, nrows) {
    .Call(`_rem_createNullEvents`, eventID, sender, target, eventAttribute, time, start, end, allEventTimes, nrows)
}

absoluteDiffAverageWeightEventAttributeCpp <- function(sender, target, time, weightvar, eventattributevar, eventattribute, xlog) {
    .Call(`_rem_absoluteDiffAverageWeightEventAttributeCpp`, sender, target, time, weightvar, eventattributevar, eventattribute, xlog)
}

