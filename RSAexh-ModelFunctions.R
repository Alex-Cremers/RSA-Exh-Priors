###############
# Vanilla RSA
###############

# S1 speakers (log scale)
RSA_S1_a_wa = function(p,lambda,deltaAB,deltaAnB){
  -log1p(exp(-lambda*(log(1-p)+deltaAnB)))
}
RSA_S1_anb_wa = function(p,lambda,deltaAB,deltaAnB){
  -log1p(exp(lambda*(log(1-p)+deltaAnB)))
}
RSA_S1_a_wab = function(p,lambda,deltaAB,deltaAnB){
  -log1p(exp(-lambda*(log(p)+deltaAB)))
}
RSA_S1_ab_wab = function(p,lambda,deltaAB,deltaAnB){
  -log1p(exp(lambda*(log(p)+deltaAB)))
}

# L1 listener
RSA_L1_a <- function(p,lambda,deltaAB,deltaAnB){
  1/(1+(1-p)*(1+exp(-lambda*(log(p)+deltaAB)))/(p*(1+exp(-lambda*(log(1-p)+deltaAnB)))))
}
RSA_L1_ab <- function(p,...){rep(1,length(p))}
RSA_L1_anb <- function(p,...){rep(0,length(p))}

# S2 speaker
RSA_S2_a_wa = function(p,lambda,deltaAB,deltaAnB){
  -log1p(exp(-lambda*(log(1-RSA_L1_a(p,lambda,deltaAB,deltaAnB))+deltaAnB)))
}
RSA_S2_anb_wa = function(p,lambda,deltaAB,deltaAnB){
  -log1p(exp(lambda*(log(1-RSA_L1_a(p,lambda,deltaAB,deltaAnB))+deltaAnB)))
}
RSA_S2_a_wab = function(p,lambda,deltaAB,deltaAnB){
  -log1p(exp(-lambda*(log(RSA_L1_a(p,lambda,deltaAB,deltaAnB))+deltaAB)))
}
RSA_S2_ab_wab = function(p,lambda,deltaAB,deltaAnB){
  -log1p(exp(lambda*(log(RSA_L1_a(p,lambda,deltaAB,deltaAnB))+deltaAB)))
}


############
# Wonky RSA
############

# t1 is the real prior
# t2 is the wonky prior
# omega = P(t=t2) (wonkiness prior)

# wRSA first-level speakers (just a prerequisite for the rest)
wRSA_S1_a_wa_t1 = function(p,lambda,deltaAB,deltaAnB){
  1/(1+exp(-lambda*(log(1-p)+deltaAnB)))
}
wRSA_S1_anb_wa_t1 = function(p,lambda,deltaAB,deltaAnB){
  1/(1+exp(lambda*(log(1-p)+deltaAnB)))
}
wRSA_S1_a_wab_t1 = function(p,lambda,deltaAB,deltaAnB){
  1/(1+exp(-lambda*(log(p)+deltaAB)))
}
wRSA_S1_ab_wab_t1 = function(p,lambda,deltaAB,deltaAnB){
  1/(1+exp(lambda*(log(p)+deltaAB)))
}
wRSA_S1_a_wa_t2 = function(p,lambda,deltaAB,deltaAnB){
  1/(1+exp(lambda*(log(2)-deltaAnB)))
}
wRSA_S1_anb_wa_t2 = function(p,lambda,deltaAB,deltaAnB){
  1/(1+exp(-lambda*(log(2)-deltaAnB)))
}
wRSA_S1_a_wab_t2 = function(p,lambda,deltaAB,deltaAnB){
  1/(1+exp(lambda*(log(2)-deltaAB)))
}
wRSA_S1_ab_wab_t2 = function(p,lambda,deltaAB,deltaAnB){
  1/(1+exp(-lambda*(log(2)-deltaAB)))
}

# wRSA listener (for ab, it's just 1):
wRSA_L1_a = function(p,lambda,deltaAB,deltaAnB,omega){
  pmin( # Use pmin to ensure we never go over 1 with rounding errors
    (omega*wRSA_S1_a_wab_t2(p,lambda,deltaAB,deltaAnB)/2+(1-omega)*wRSA_S1_a_wab_t1(p,lambda,deltaAB,deltaAnB)*p)/
    (omega*wRSA_S1_a_wab_t2(p,lambda,deltaAB,deltaAnB)/2+(1-omega)*p*wRSA_S1_a_wab_t1(p,lambda,deltaAB,deltaAnB)+
       omega*wRSA_S1_a_wa_t2(p,lambda,deltaAB,deltaAnB)/2+(1-omega)*(1-p)*wRSA_S1_a_wa_t1(p,lambda,deltaAB,deltaAnB)),
    1)
}

# wRSA posterior on wonkiness, only here for better understanding of the model.
wRSA_L1_t_a = function(p,lambda,deltaAB,deltaAnB,omega){
  omega*0.5*(wRSA_S1_a_wab_t2(p,lambda,deltaAB,deltaAnB)+wRSA_S1_a_wa_t2(p,lambda,deltaAB,deltaAnB))/
    (omega*wRSA_S1_a_wab_t2(p,lambda,deltaAB,deltaAnB)/2+(1-omega)*p*wRSA_S1_a_wab_t1(p,lambda,deltaAB,deltaAnB)+
       omega*wRSA_S1_a_wa_t2(p,lambda,deltaAB,deltaAnB)/2+(1-omega)*(1-p)*wRSA_S1_a_wa_t1(p,lambda,deltaAB,deltaAnB))
}


# wRSA speakers on log scale (use log1p for precision)
wRSA_S2_a_wa = function(p,lambda,deltaAB,deltaAnB,omega){
  -log1p(exp(-lambda*(log(1-wRSA_L1_a(p,lambda,deltaAB,deltaAnB,omega)) + deltaAnB)))
}
wRSA_S2_anb_wa = function(p,lambda,deltaAB,deltaAnB,omega){
  -log1p(exp(lambda*(log(1-wRSA_L1_a(p,lambda,deltaAB,deltaAnB,omega)) + deltaAnB)))
}
wRSA_S2_a_wab = function(p,lambda,deltaAB,deltaAnB,omega){
  -log1p(exp(-lambda*(log(wRSA_L1_a(p,lambda,deltaAB,deltaAnB,omega)) + deltaAB)))
}
wRSA_S2_ab_wab = function(p,lambda,deltaAB,deltaAnB,omega){
  -log1p(exp(lambda*(log(wRSA_L1_a(p,lambda,deltaAB,deltaAnB,omega)) + deltaAB)))
}


#################
# Bayesian wonky
#################


# BwRSA listener (for ab, it's just 1):
BwRSA_L1_a = function(p,lambda,deltaAB,deltaAnB,omega){
  pmin( # use pmin to ensure that rounding errors don't push us above 1
    p*(omega*wRSA_S1_a_wab_t2(p,lambda,deltaAB,deltaAnB)+(1-omega)*wRSA_S1_a_wab_t1(p,lambda,deltaAB,deltaAnB))/
    (omega*p*wRSA_S1_a_wab_t2(p,lambda,deltaAB,deltaAnB)+(1-omega)*p*wRSA_S1_a_wab_t1(p,lambda,deltaAB,deltaAnB)+
       omega*(1-p)*wRSA_S1_a_wa_t2(p,lambda,deltaAB,deltaAnB)+(1-omega)*(1-p)*wRSA_S1_a_wa_t1(p,lambda,deltaAB,deltaAnB)),
    1)
}


# BwRSA speakers on log scale (use log1p for precision)
BwRSA_S2_a_wa = function(p,lambda,deltaAB,deltaAnB,omega){
  -log1p(exp(-lambda*(log(1-BwRSA_L1_a(p,lambda,deltaAB,deltaAnB,omega)) + deltaAnB)))
}
BwRSA_S2_anb_wa = function(p,lambda,deltaAB,deltaAnB,omega){
  -log1p(exp(lambda*(log(1-BwRSA_L1_a(p,lambda,deltaAB,deltaAnB,omega)) + deltaAnB)))
}
BwRSA_S2_a_wab = function(p,lambda,deltaAB,deltaAnB,omega){
  -log1p(exp(-lambda*(log(BwRSA_L1_a(p,lambda,deltaAB,deltaAnB,omega)) + deltaAB)))
}
BwRSA_S2_ab_wab = function(p,lambda,deltaAB,deltaAnB,omega){
  -log1p(exp(lambda*(log(BwRSA_L1_a(p,lambda,deltaAB,deltaAnB,omega)) + deltaAB)))
}

#################
#   svRSA model 
# Spector (2017)
#################

# Qa is the partial QUD
# Qb is the total QUD
# q = P(Qb) (prior on QUD)
# chi = P(exh) (prior on interpretation, to be set to 1/2)

svRSA_S1_a_Qa <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  1/(1+exp(-lambda*deltaAnB)+exp(-lambda*deltaAB))
}
svRSA_S1_anb_Qa <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  1/(1+exp(lambda*deltaAnB)+exp(lambda*(deltaAnB-deltaAB)))
}
svRSA_S1_ab_Qa <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  1/(1+exp(lambda*deltaAB)+exp(lambda*(deltaAB-deltaAnB)))
}
svRSA_S1_a_wa_Qb <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  1/(1+exp(-lambda*((1-chi)*log(1-p)+deltaAnB)))
}
svRSA_S1_anb_wa_Qb <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  1/(1+exp(lambda*((1-chi)*log(1-p)+deltaAnB)))
}

# svRSA L1 functions (full)
svRSA_L1_Qa_wa_a <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  (1-p)*(1-q)*svRSA_S1_a_Qa(p,lambda,deltaAB,deltaAnB,q,chi)/
    ((1-q)*svRSA_S1_a_Qa(p,lambda,deltaAB,deltaAnB,q,chi) + q*(1-p)*svRSA_S1_a_wa_Qb(p,lambda,deltaAB,deltaAnB,q,chi))
}
svRSA_L1_Qa_wab_a <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  p*(1-q)*svRSA_S1_a_Qa(p,lambda,deltaAB,deltaAnB,q,chi)/
    ((1-q)*svRSA_S1_a_Qa(p,lambda,deltaAB,deltaAnB,q,chi) + q*(1-p)*svRSA_S1_a_wa_Qb(p,lambda,deltaAB,deltaAnB,q,chi))
}
svRSA_L1_Qa_wa_anb <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  (1-p)*(1-q)*svRSA_S1_anb_Qa(p,lambda,deltaAB,deltaAnB,q,chi)/
    ((1-q)*svRSA_S1_anb_Qa(p,lambda,deltaAB,deltaAnB,q,chi) + q*(1-p)*svRSA_S1_anb_wa_Qb(p,lambda,deltaAB,deltaAnB,q,chi))
}
svRSA_L1_Qa_wab_anb <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  p*(1-q)*svRSA_S1_anb_Qa(p,lambda,deltaAB,deltaAnB,q,chi)/
    ((1-q)*svRSA_S1_anb_Qa(p,lambda,deltaAB,deltaAnB,q,chi) + q*(1-p)*svRSA_S1_anb_wa_Qb(p,lambda,deltaAB,deltaAnB,q,chi))
}
svRSA_L1_Qa_wa_ab <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  (1-p)*(1-q)*svRSA_S1_ab_Qa(p,lambda,deltaAB,deltaAnB,q,chi)/
    ((1-q)*svRSA_S1_ab_Qa(p,lambda,deltaAB,deltaAnB,q,chi) + q*p)
}
svRSA_L1_Qa_wab_ab <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  p*(1-q)*svRSA_S1_ab_Qa(p,lambda,deltaAB,deltaAnB,q,chi)/
    ((1-q)*svRSA_S1_ab_Qa(p,lambda,deltaAB,deltaAnB,q,chi) + q*p)
}
svRSA_L1_Qb_wa_a <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  (1-p)*q*svRSA_S1_a_wa_Qb(p,lambda,deltaAB,deltaAnB,q,chi)/
    ((1-q)*svRSA_S1_a_Qa(p,lambda,deltaAB,deltaAnB,q,chi) + q*(1-p)*svRSA_S1_a_wa_Qb(p,lambda,deltaAB,deltaAnB,q,chi))
}
svRSA_L1_Qb_wa_anb <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  (1-p)*q*svRSA_S1_anb_wa_Qb(p,lambda,deltaAB,deltaAnB,q,chi)/
    ((1-q)*svRSA_S1_anb_Qa(p,lambda,deltaAB,deltaAnB,q,chi) + q*(1-p)*svRSA_S1_anb_wa_Qb(p,lambda,deltaAB,deltaAnB,q,chi))
}
svRSA_L1_Qb_wab_ab <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  p*q/((1-q)*svRSA_S1_ab_Qa(p,lambda,deltaAB,deltaAnB,q,chi) + q*p)
}


# svRSA logL1 functions (marginalized Qa posterior, useful for the speaker who conveys Qa, since she does not care about w)
svRSA_L1_Qa_a <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  log((1-q)*svRSA_S1_a_Qa(p,lambda,deltaAB,deltaAnB,q,chi)/((1-q)*svRSA_S1_a_Qa(p,lambda,deltaAB,deltaAnB,q,chi)+q*(1-p)*svRSA_S1_a_wa_Qb(p,lambda,deltaAB,deltaAnB,q,chi)))
}
svRSA_L1_Qa_ab <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  log((1-q)*svRSA_S1_ab_Qa(p,lambda,deltaAB,deltaAnB,q,chi)/(q*p+(1-q)*svRSA_S1_ab_Qa(p,lambda,deltaAB,deltaAnB,q,chi)))
}
svRSA_L1_Qa_anb <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  log((1-q)*svRSA_S1_anb_Qa(p,lambda,deltaAB,deltaAnB,q,chi)/((1-q)*svRSA_S1_anb_Qa(p,lambda,deltaAB,deltaAnB,q,chi)+q*(1-p)*svRSA_S1_anb_wa_Qb(p,lambda,deltaAB,deltaAnB,q,chi)))
}
# svRSA diff logL1 functions (to simplify S2 computations)
svRSA_L1_Qa_aab <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  log(svRSA_S1_a_Qa(p,lambda,deltaAB,deltaAnB,q,chi)*
        (q*p+(1-q)*svRSA_S1_ab_Qa(p,lambda,deltaAB,deltaAnB,q,chi))/
        (svRSA_S1_ab_Qa(p,lambda,deltaAB,deltaAnB,q,chi)*
           ((1-q)*svRSA_S1_a_Qa(p,lambda,deltaAB,deltaAnB,q,chi)+
              q*(1-p)*svRSA_S1_a_wa_Qb(p,lambda,deltaAB,deltaAnB,q,chi)
            )
         )
      )
}
svRSA_L1_Qa_aanb <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  log(svRSA_S1_a_Qa(p,lambda,deltaAB,deltaAnB,q,chi)*
        ((1-q)*svRSA_S1_anb_Qa(p,lambda,deltaAB,deltaAnB,q,chi)+q*(1-p)*svRSA_S1_anb_wa_Qb(p,lambda,deltaAB,deltaAnB,q,chi))/
        (svRSA_S1_anb_Qa(p,lambda,deltaAB,deltaAnB,q,chi)*
           ((1-q)*svRSA_S1_a_Qa(p,lambda,deltaAB,deltaAnB,q,chi)+
              q*(1-p)*svRSA_S1_a_wa_Qb(p,lambda,deltaAB,deltaAnB,q,chi)
            )
         )
        )
}
svRSA_L1_Qa_abanb <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  log(svRSA_S1_ab_Qa(p,lambda,deltaAB,deltaAnB,q,chi)*
        ((1-q)*svRSA_S1_anb_Qa(p,lambda,deltaAB,deltaAnB,q,chi)+q*(1-p)*svRSA_S1_anb_wa_Qb(p,lambda,deltaAB,deltaAnB,q,chi))/
        (svRSA_S1_anb_Qa(p,lambda,deltaAB,deltaAnB,q,chi)*
           (q*p+(1-q)*svRSA_S1_ab_Qa(p,lambda,deltaAB,deltaAnB,q,chi))
         )
        )
}

# svRSA L1 functions (marginalized worlds posterior for comprehension data)
svRSA_L1_a <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  svRSA_L1_Qa_wab_a(p,lambda,deltaAB,deltaAnB,q,chi)
}
svRSA_L1_ab <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  svRSA_L1_Qb_wab_ab(p,lambda,deltaAB,deltaAnB,q,chi)+
    svRSA_L1_Qa_wab_ab(p,lambda,deltaAB,deltaAnB,q,chi)
}
svRSA_L1_anb <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  svRSA_L1_Qa_wab_anb(p,lambda,deltaAB,deltaAnB,q,chi)
}

# svRSA log-S2 functions, conveying world and QUD Qb
log_svRSA_S2_a_wa_Qb <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  ifelse(p==1,-Inf,-log1p(
    exp(lambda*(log(svRSA_L1_Qb_wa_anb(p,lambda,deltaAB,deltaAnB,q,chi))-log(svRSA_L1_Qb_wa_a(p,lambda,deltaAB,deltaAnB,q,chi))-deltaAnB))))
}
log_svRSA_S2_anb_wa_Qb <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  ifelse(p==1,0,-log1p(
    exp(-lambda*(log(svRSA_L1_Qb_wa_anb(p,lambda,deltaAB,deltaAnB,q,chi))-log(svRSA_L1_Qb_wa_a(p,lambda,deltaAB,deltaAnB,q,chi))-deltaAnB))))
}


# svRSA log-S2 functions, conveying (world and) QUD Qa
log_svRSA_S2_a_Qa <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  -log1p(
    exp(-lambda*(svRSA_L1_Qa_aanb(p,lambda,deltaAB,deltaAnB,q,chi)+deltaAnB))+
    exp(-lambda*(svRSA_L1_Qa_aab(p,lambda,deltaAB,deltaAnB,q,chi)+deltaAB))
    )
}
log_svRSA_S2_ab_Qa <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  -log1p(
    exp(lambda*(-svRSA_L1_Qa_abanb(p,lambda,deltaAB,deltaAnB,q,chi)+deltaAB-deltaAnB))+
      exp(lambda*(svRSA_L1_Qa_aab(p,lambda,deltaAB,deltaAnB,q,chi)+deltaAB))
  )
}
log_svRSA_S2_anb_Qa <- function(p,lambda,deltaAB,deltaAnB,q,chi){
  -log1p(
    exp(lambda*(svRSA_L1_Qa_abanb(p,lambda,deltaAB,deltaAnB,q,chi)-deltaAB+deltaAnB))+
      exp(lambda*(svRSA_L1_Qa_aanb(p,lambda,deltaAB,deltaAnB,q,chi)+deltaAnB))
  )
}


# Mixture S2 for svRSA1 model S2 : P(Qa)S2(|Qa)+P(Qb)S2(|Qb) with chi=1/2
svRSA1_S2_ab_wab <- function(p,lambda,deltaAB,deltaAnB,q){
  log(q+(1-q)*exp(log_svRSA_S2_ab_Qa(p,lambda,deltaAB,deltaAnB,q,0.5)))
}
svRSA1_S2_ab_wa <- function(p,lambda,deltaAB,deltaAnB,q){
  log(1-q)+log_svRSA_S2_ab_Qa(p,lambda,deltaAB,deltaAnB,q,0.5)
}
svRSA1_S2_a_wa <- function(p,lambda,deltaAB,deltaAnB,q){
  log(q*exp(log_svRSA_S2_a_wa_Qb(p,lambda,deltaAB,deltaAnB,q,0.5))+(1-q)*exp(log_svRSA_S2_a_Qa(p,lambda,deltaAB,deltaAnB,q,0.5)))
}
svRSA1_S2_a_wab <- function(p,lambda,deltaAB,deltaAnB,q){
  log(1-q) + log_svRSA_S2_a_Qa(p,lambda,deltaAB,deltaAnB,q,0.5)
}
svRSA1_S2_anb_wab <- function(p,lambda,deltaAB,deltaAnB,q){
  log(1-q) + log_svRSA_S2_anb_Qa(p,lambda,deltaAB,deltaAnB,q,0.5)
}
svRSA1_S2_anb_wa <- function(p,lambda,deltaAB,deltaAnB,q){
  log(q*exp(log_svRSA_S2_anb_wa_Qb(p,lambda,deltaAB,deltaAnB,q,0.5))+(1-q)*exp(log_svRSA_S2_anb_Qa(p,lambda,deltaAB,deltaAnB,q,0.5)))
}

######################
# Lexical uncertainty
######################

# Possible reinforcements for 'a':
# * a_1 = literal
# * a_2 = exhaustive
# * a_3 = anti-exhaustive

# We start with a generic LU model with free priors on lexica:
# lex_prior = c(P(a1),P(a2),P(a3)) (must sum to 1)
# We can get either FREE-LU or EXH-LU by setting the priors differently

RSA_LU_S1_wa_a_1 <- function(p,lambda,deltaAB,deltaAnB){1/(1+exp(-lambda*(log(1-p)+deltaAnB)))}
RSA_LU_S1_wa_a_2 <- function(p,lambda,deltaAB,deltaAnB){1/(1+exp(-lambda*deltaAnB))}
RSA_LU_S1_wab_a_1 <- function(p,lambda,deltaAB,deltaAnB){1/(1+exp(-lambda*(log(p)+deltaAB)))}
RSA_LU_S1_wab_a_3 <- function(p,lambda,deltaAB,deltaAnB){1/(1+exp(-lambda*deltaAB))}

RSA_LU_L1_a <- function(p,lambda,deltaAB,deltaAnB,lex_prior){
  p*(lex_prior[1]*RSA_LU_S1_wab_a_1(p,lambda,deltaAB,deltaAnB)+lex_prior[3]*RSA_LU_S1_wab_a_3(p,lambda,deltaAB,deltaAnB))/
    (p*lex_prior[1]*RSA_LU_S1_wab_a_1(p,lambda,deltaAB,deltaAnB)+p*lex_prior[3]*RSA_LU_S1_wab_a_3(p,lambda,deltaAB,deltaAnB)+
       (1-p)*lex_prior[1]*RSA_LU_S1_wa_a_1(p,lambda,deltaAB,deltaAnB)+(1-p)*lex_prior[2]*RSA_LU_S1_wa_a_2(p,lambda,deltaAB,deltaAnB))
}
RSA_LU_S2_a_wab <- function(p,lambda,deltaAB,deltaAnB,lex_prior){
  -log1p(exp(-lambda*(log(RSA_LU_L1_a(p,lambda,deltaAB,deltaAnB,lex_prior))+deltaAB)))
}
RSA_LU_S2_ab_wab <- function(p,lambda,deltaAB,deltaAnB,lex_prior){
  -log1p(exp(lambda*(log(RSA_LU_L1_a(p,lambda,deltaAB,deltaAnB,lex_prior))+deltaAB)))
}
RSA_LU_S2_a_wa <- function(p,lambda,deltaAB,deltaAnB,lex_prior){
  -log1p(exp(-lambda*(log(1-RSA_LU_L1_a(p,lambda,deltaAB,deltaAnB,lex_prior))+deltaAnB)))
}
RSA_LU_S2_anb_wa <- function(p,lambda,deltaAB,deltaAnB,lex_prior){
  -log1p(exp(lambda*(log(1-RSA_LU_L1_a(p,lambda,deltaAB,deltaAnB,lex_prior))+deltaAnB)))
}


#######################
#   Lexical Intensions
# Franke&Bergen (2019)
#######################

# Marginal S1 speaker
RSA_LI_S1_a_wa<- function(p,lambda,deltaAB,deltaAnB){
  log1p(exp(lambda*log(1-p)))-log1p(exp(lambda*log(1-p))+2*exp(-lambda*deltaAnB))
}
RSA_LI_S1_a_wab<- function(p,lambda,deltaAB,deltaAnB){
  -log1p(2*exp(-lambda*(log(p)+deltaAB)))
}
RSA_LI_S1_ab_wab<- function(p,lambda,deltaAB,deltaAnB){
  -log1p(exp(lambda*log(p)+lambda*deltaAB-log(2)))
}
RSA_LI_S1_anb_wa<- function(p,lambda,deltaAB,deltaAnB){
  -log1p(exp(lambda*deltaAnB-log(2))+exp(lambda*(log(1-p)+deltaAnB)-log(2)))
}

# Marginal L1 listener, logL1(w_ab|a) and logL1(w_a|a) (for use in S2)
RSA_LI_L1_a <- function(p,lambda,deltaAB,deltaAnB){
  1/(1+(1-p)*(1+2*exp(-lambda*(log(p)+deltaAB)))*(1+exp(lambda*log(1-p)))/(p*(1+exp(lambda*log(1-p))+2*exp(-lambda*deltaAnB))))
}
RSA_LI_logL1_a <- function(p,lambda,deltaAB,deltaAnB){
  -log1p((1-p)*(1+2*exp(-lambda*(log(p)+deltaAB)))*(1+exp(lambda*log(1-p)))/(p*(1+exp(lambda*log(1-p))+2*exp(-lambda*deltaAnB))))
}
RSA_LI_logcL1_a <- function(p,lambda,deltaAB,deltaAnB){
  -log1p(p*(1+exp(lambda*log(1-p))+2*exp(-lambda*deltaAnB))/((1-p)*(1+2*exp(-lambda*(log(p)+deltaAB)))*(1+exp(lambda*log(1-p)))))
}

# S2 speaker
RSA_LI_S2_a_wa <- function(p,lambda,deltaAB,deltaAnB){
  -log1p(exp(-lambda*(RSA_LI_logcL1_a(p,lambda,deltaAB,deltaAnB)+deltaAnB)))
}
RSA_LI_S2_anb_wa <- function(p,lambda,deltaAB,deltaAnB){
  -log1p(exp(lambda*(RSA_LI_logcL1_a(p,lambda,deltaAB,deltaAnB)+deltaAnB)))
}
RSA_LI_S2_a_wab <- function(p,lambda,deltaAB,deltaAnB){
  -log1p(exp(-lambda*(RSA_LI_logL1_a(p,lambda,deltaAB,deltaAnB)+deltaAB)))
}
RSA_LI_S2_ab_wab <- function(p,lambda,deltaAB,deltaAnB){
  -log1p(exp(lambda*(RSA_LI_logL1_a(p,lambda,deltaAB,deltaAnB)+deltaAB)))
}

#################################################
# Define the likelihood functions for each model
#################################################

RSA_Production_LogLikelihood <- function(Message,World,Prior,lambda,deltaAB,deltaAnB){
  case_when(
    Message == "a&b" & World == "ab" ~ RSA_S2_ab_wab(Prior,lambda,deltaAB,deltaAnB),
    Message %in% c("a&!b","!b") & World == "a" ~ RSA_S2_anb_wa(Prior,lambda,deltaAB,deltaAnB),
    Message  == "a" & World == "a" ~ RSA_S2_a_wa(Prior,lambda,deltaAB,deltaAnB),
    Message  == "a" & World == "ab" ~ RSA_S2_a_wab(Prior,lambda,deltaAB,deltaAnB),
    T ~ -Inf
  )
}
RSA_Comprehension_LogLikelihood <- function(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB,sigma_a,sigma_ab){
  L1a=RSA_L1_a(Prior,lambda,deltaAB,deltaAnB)
  case_when(
    Posterior==1 ~ pnorm(1,if_else(Utterance=="a",L1a,1),if_else(Utterance=="a",sigma_a,sigma_ab),lower.tail = F,log.p = T),
    Posterior==0 ~ pnorm(0,if_else(Utterance=="a",L1a,1),if_else(Utterance=="a",sigma_a,sigma_ab),lower.tail = T,log.p = T),
    T ~ dnorm(Posterior,if_else(Utterance=="a",L1a,1),if_else(Utterance=="a",sigma_a,sigma_ab),log = T)
  )
}

wRSA_Production_LogLikelihood <- function(Message,World,Prior,lambda,deltaAB,deltaAnB,omega){
  case_when(
    Message == "a&b" & World == "ab" ~ wRSA_S2_ab_wab(Prior,lambda,deltaAB,deltaAnB,omega),
    Message %in% c("a&!b","!b") & World == "a" ~ wRSA_S2_anb_wa(Prior,lambda,deltaAB,deltaAnB,omega),
    Message  == "a" & World == "a" ~ wRSA_S2_a_wa(Prior,lambda,deltaAB,deltaAnB,omega),
    Message  == "a" & World == "ab" ~ wRSA_S2_a_wab(Prior,lambda,deltaAB,deltaAnB,omega),
    T ~ -Inf
  )
}
wRSA_Comprehension_LogLikelihood <- function(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB,omega,sigma_a,sigma_ab){
  L1a=wRSA_L1_a(Prior,lambda,deltaAB,deltaAnB,omega)
  case_when(
    Posterior==1 ~ pnorm(1,if_else(Utterance=="a",L1a,1),if_else(Utterance=="a",sigma_a,sigma_ab),lower.tail = F,log.p = T),
    Posterior==0 ~ pnorm(0,if_else(Utterance=="a",L1a,1),if_else(Utterance=="a",sigma_a,sigma_ab),lower.tail = T,log.p = T),
    T ~ dnorm(Posterior,if_else(Utterance=="a",L1a,1),if_else(Utterance=="a",sigma_a,sigma_ab),log = T)
  )
}

BwRSA_Production_LogLikelihood <- function(Message,World,Prior,lambda,deltaAB,deltaAnB,omega){
  case_when(
    Message == "a&b" & World == "ab" ~ BwRSA_S2_ab_wab(Prior,lambda,deltaAB,deltaAnB,omega),
    Message %in% c("a&!b","!b") & World == "a" ~ BwRSA_S2_anb_wa(Prior,lambda,deltaAB,deltaAnB,omega),
    Message  == "a" & World == "a" ~ BwRSA_S2_a_wa(Prior,lambda,deltaAB,deltaAnB,omega),
    Message  == "a" & World == "ab" ~ BwRSA_S2_a_wab(Prior,lambda,deltaAB,deltaAnB,omega),
    T ~ -Inf
  )
}
BwRSA_Comprehension_LogLikelihood <- function(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB,omega,sigma_a,sigma_ab){
  L1a=BwRSA_L1_a(Prior,lambda,deltaAB,deltaAnB,omega)
  case_when(
    Posterior==1 ~ pnorm(1,if_else(Utterance=="a",L1a,1),if_else(Utterance=="a",sigma_a,sigma_ab),lower.tail = F,log.p = T),
    Posterior==0 ~ pnorm(0,if_else(Utterance=="a",L1a,1),if_else(Utterance=="a",sigma_a,sigma_ab),lower.tail = T,log.p = T),
    T ~ dnorm(Posterior,if_else(Utterance=="a",L1a,1),if_else(Utterance=="a",sigma_a,sigma_ab),log = T)
  )
}

# svRSA with P(Qa)S2(|Qa)+P(Qb)S2(|Qb)
svRSA1_Production_LogLikelihood <- function(Message,World,Prior,lambda,deltaAB,deltaAnB,q){
  case_when(
    Message == "a&b" & World == "ab"~ svRSA1_S2_ab_wab(Prior,lambda,deltaAB,deltaAnB,q),
    Message == "a&b" & World == "a"~ svRSA1_S2_a_wab(Prior,lambda,deltaAB,deltaAnB,q),
    Message %in% c("a&!b","!b") & World == "a"~ svRSA1_S2_anb_wa(Prior,lambda,deltaAB,deltaAnB,q),
    Message %in% c("a&!b","!b") & World == "ab"~ svRSA1_S2_anb_wab(Prior,lambda,deltaAB,deltaAnB,q),
    Message  == "a" & World == "a" ~ svRSA1_S2_a_wa(Prior,lambda,deltaAB,deltaAnB,q),
    Message  == "a" & World == "ab" ~ svRSA1_S2_a_wab(Prior,lambda,deltaAB,deltaAnB,q),
    T ~ -Inf
  )
}

# svRSA with S2(|Qb)
svRSA2_Production_LogLikelihood <- function(Message,World,Prior,lambda,deltaAB,deltaAnB,q){
  case_when(
    Message == "a&b" & World == "a"~ -Inf,
    Message == "a&b" & World == "ab"~ 0,
    Message %in% c("a&!b","!b") & World == "a"~ log_svRSA_S2_anb_wa_Qb(Prior,lambda,deltaAB,deltaAnB,q,0.5),
    Message %in% c("a&!b","!b") & World == "ab"~ -Inf,
    Message  == "a" & World == "a" ~ log_svRSA_S2_a_wa_Qb(Prior,lambda,deltaAB,deltaAnB,q,0.5),
    Message  == "a" & World == "ab" ~ -Inf,
    T ~ -Inf
  )
}


svRSA_Comprehension_LogLikelihood <- function(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB,q,sigma_a,sigma_ab){
  L1a=svRSA_L1_a(Prior,lambda,deltaAB,deltaAnB,q,0.5)
  L1ab=svRSA_L1_ab(Prior,lambda,deltaAB,deltaAnB,q,0.5)
  case_when(
    Posterior==1 ~ pnorm(1,if_else(Utterance=="a",L1a,L1ab),if_else(Utterance=="a",sigma_a,sigma_ab),lower.tail = F,log.p = T),
    Posterior==0 ~ pnorm(0,if_else(Utterance=="a",L1a,L1ab),if_else(Utterance=="a",sigma_a,sigma_ab),lower.tail = T,log.p = T),
    T ~ dnorm(Posterior,if_else(Utterance=="a",L1a,L1ab),if_else(Utterance=="a",sigma_a,sigma_ab),log = T)
  )
}
# For compatibility:
svRSA1_Comprehension_LogLikelihood <- function(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB,q,sigma_a,sigma_ab){svRSA_Comprehension_LogLikelihood(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB,q,sigma_a,sigma_ab)}
svRSA2_Comprehension_LogLikelihood <- function(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB,q,sigma_a,sigma_ab){svRSA_Comprehension_LogLikelihood(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB,q,sigma_a,sigma_ab)}

# RSA-LU functions parametrized as follows:
# p1 is the probability of lexical strengthening
# p2 is the probability of "a" being interpreted as "a&!b" given strengthening.
# Uniform prior: p1=2/3, p2=1/2
# Grammatical prior: p1=1/2, p2=1

RSA_LU_Production_LogLikelihood <- function(Message,World,Prior,lambda,deltaAB,deltaAnB,p1,p2){
  lex_prior <- c(1-p1,p1*p2,p1*(1-p2)) # literal, exhaustive, anti-exhaustive
  case_when(
    Message == "a&b" & World == "ab" ~ RSA_LU_S2_ab_wab(Prior,lambda,deltaAB,deltaAnB,lex_prior),
    Message %in% c("a&!b","!b") & World == "a" ~ RSA_LU_S2_anb_wa(Prior,lambda,deltaAB,deltaAnB,lex_prior),
    Message  == "a" & World == "a" ~ RSA_LU_S2_a_wa(Prior,lambda,deltaAB,deltaAnB,lex_prior),
    Message  == "a" & World == "ab" ~ RSA_LU_S2_a_wab(Prior,lambda,deltaAB,deltaAnB,lex_prior),
    T ~ -Inf
  )
}
RSA_LU_Comprehension_LogLikelihood <- function(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB,p1,p2,sigma_a,sigma_ab){
  lex_prior <- c(1-p1,p1*p2,p1*(1-p2))
  L1a=RSA_LU_L1_a(Prior,lambda,deltaAB,deltaAnB,lex_prior)
  case_when(
    Posterior==1 ~ pnorm(1,if_else(Utterance=="a",L1a,1),if_else(Utterance=="a",sigma_a,sigma_ab),lower.tail = F,log.p = T),
    Posterior==0 ~ pnorm(0,if_else(Utterance=="a",L1a,1),if_else(Utterance=="a",sigma_a,sigma_ab),lower.tail = T,log.p = T),
    T ~ dnorm(Posterior,if_else(Utterance=="a",L1a,1),if_else(Utterance=="a",sigma_a,sigma_ab),log = T)
  )
}


# FREE_LU : uniform priors on all possible reinforcements:
FREE_LU_Comprehension_LogLikelihood <- function(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB,sigma_a,sigma_ab){
  RSA_LU_Comprehension_LogLikelihood(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB,2/3,1/2,sigma_a,sigma_ab)
}
FREE_LU_Production_LogLikelihood <- function(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB){
  RSA_LU_Production_LogLikelihood(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB,2/3,1/2)
}
# EXH_LU : uniform priors but symmetry broken in the grammar:
EXH_LU_Comprehension_LogLikelihood <- function(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB,sigma_a,sigma_ab){
  RSA_LU_Comprehension_LogLikelihood(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB,1/2,1,sigma_a,sigma_ab)
}
EXH_LU_Production_LogLikelihood <- function(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB){
  RSA_LU_Production_LogLikelihood(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB,1/2,1)
}


# LI model using S2 speaker
RSA_LI2_S2_a_wa <- RSA_LI_S2_a_wa
RSA_LI2_S2_a_wab <- RSA_LI_S2_a_wab
RSA_LI2_S2_ab_wab <- RSA_LI_S2_ab_wab
RSA_LI2_S2_anb_wa <- RSA_LI_S2_anb_wa
RSA_LI2_Production_LogLikelihood <- function(Message,World,Prior,lambda,deltaAB,deltaAnB){
  case_when(
    Message == "a&b" & World == "ab" ~ RSA_LI_S2_ab_wab(Prior,lambda,deltaAB,deltaAnB),
    Message %in% c("a&!b","!b") & World == "a" ~ RSA_LI_S2_anb_wa(Prior,lambda,deltaAB,deltaAnB),
    Message  == "a" & World == "a" ~ RSA_LI_S2_a_wa(Prior,lambda,deltaAB,deltaAnB),
    Message  == "a" & World == "ab" ~ RSA_LI_S2_a_wab(Prior,lambda,deltaAB,deltaAnB),
    T ~ -Inf
  )
}

RSA_LI2_Comprehension_LogLikelihood <- function(Posterior,Utterance,Prior,lambda,deltaAB,deltaAnB,sigma_a,sigma_ab){
  L1a=RSA_LI_L1_a(Prior,lambda,deltaAB,deltaAnB)
  case_when(
    Posterior==1 ~ pnorm(1,if_else(Utterance=="a",L1a,1),if_else(Utterance=="a",sigma_a,sigma_ab),lower.tail = F,log.p = T),
    Posterior==0 ~ pnorm(0,if_else(Utterance=="a",L1a,1),if_else(Utterance=="a",sigma_a,sigma_ab),lower.tail = T,log.p = T),
    T ~ dnorm(Posterior,if_else(Utterance=="a",L1a,1),if_else(Utterance=="a",sigma_a,sigma_ab),log = T)
  )
}

# LI model using marginal S1 speaker
# (must still name it S2 for compatibility with main script functions)
RSA_LI1_S2_a_wa <- RSA_LI_S1_a_wa
RSA_LI1_S2_a_wab <- RSA_LI_S1_a_wab
RSA_LI1_S2_ab_wab <- RSA_LI_S1_ab_wab
RSA_LI1_S2_anb_wa <- RSA_LI_S1_anb_wa
RSA_LI1_Production_LogLikelihood <- function(Message,World,Prior,lambda,deltaAB,deltaAnB){
  case_when(
    Message == "a&b" & World == "ab" ~ RSA_LI1_S2_ab_wab(Prior,lambda,deltaAB,deltaAnB),
    Message %in% c("a&!b","!b") & World == "a" ~ RSA_LI1_S2_anb_wa(Prior,lambda,deltaAB,deltaAnB),
    Message  == "a" & World == "a" ~ RSA_LI1_S2_a_wa(Prior,lambda,deltaAB,deltaAnB),
    Message  == "a" & World == "ab" ~ RSA_LI1_S2_a_wab(Prior,lambda,deltaAB,deltaAnB),
    T ~ -Inf
  )
}
RSA_LI1_Comprehension_LogLikelihood <- RSA_LI2_Comprehension_LogLikelihood

# Define missing listener functions for graphs:
RSA_L1_ab <- function(p,...){rep(1,length(p))}
wRSA_L1_ab <- function(p,...){rep(1,length(p))}
BwRSA_L1_ab <- function(p,...){rep(1,length(p))}
svRSA1_L1_a <- function(p,lambda,deltaAB,deltaAnB,q){svRSA_L1_a(p,lambda,deltaAB,deltaAnB,q,0.5)}
svRSA1_L1_ab <- function(p,lambda,deltaAB,deltaAnB,q){svRSA_L1_ab(p,lambda,deltaAB,deltaAnB,q,0.5)}
svRSA2_L1_a <- function(p,lambda,deltaAB,deltaAnB,q){svRSA_L1_a(p,lambda,deltaAB,deltaAnB,q,0.5)}
svRSA2_L1_ab <- function(p,lambda,deltaAB,deltaAnB,q){svRSA_L1_ab(p,lambda,deltaAB,deltaAnB,q,0.5)}
FREE_LU_L1_a <- function(p,lambda,deltaAB,deltaAnB,lex_prior){RSA_LU_L1_a(p,lambda,deltaAB,deltaAnB,rep(1/3,3))}
EXH_LU_L1_a <- function(p,lambda,deltaAB,deltaAnB,lex_prior){RSA_LU_L1_a(p,lambda,deltaAB,deltaAnB,c(1/2,1/2,0))}
FREE_LU_L1_ab <- function(p,...){rep(1,length(p))}
EXH_LU_L1_ab <- function(p,...){rep(1,length(p))}
RSA_LI1_L1_a <- RSA_LI_L1_a
RSA_LI1_L1_ab <- function(p,...){rep(1,length(p))}
RSA_LI2_L1_a <- RSA_LI_L1_a
RSA_LI2_L1_ab <- RSA_LI1_L1_ab

# Define missing speaker functions for graphs:
svRSA2_S2_a_wa <- function(p,lambda,deltaAB,deltaAnB,q){log_svRSA_S2_a_wa_Qb(p,lambda,deltaAB,deltaAnB,q,0.5)}
svRSA2_S2_ab_wab <- function(p,lambda,deltaAB,deltaAnB,q){rep(0,length(p))}
FREE_LU_S2_a_wa <- function(p,lambda,deltaAB,deltaAnB){RSA_LU_S2_a_wa(p,lambda,deltaAB,deltaAnB,rep(1/3,3))}
FREE_LU_S2_ab_wab <- function(p,lambda,deltaAB,deltaAnB){RSA_LU_S2_ab_wab(p,lambda,deltaAB,deltaAnB,rep(1/3,3))}
EXH_LU_S2_a_wa <- function(p,lambda,deltaAB,deltaAnB){RSA_LU_S2_a_wa(p,lambda,deltaAB,deltaAnB,c(1/2,1/2,0))}
EXH_LU_S2_ab_wab <- function(p,lambda,deltaAB,deltaAnB){RSA_LU_S2_ab_wab(p,lambda,deltaAB,deltaAnB,c(1/2,1/2,0))}

