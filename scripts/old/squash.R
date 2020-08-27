play_point = function(p){
  return(rbern(1, p))
}

HiHo_win = function(s1, s2){

  if (s1 == 9 & s2 < 8) return('1')
  if (s2 == 9 & s1 < 8) return('2')
  
  if (s1 == 10) return('1')
  if (s2 == 10) return('2')
  
  return(F)
}



PAR_win = function(s1, s2){
  
  if (s1 >= 11 & s1-s2 >=2) return('1')
  if (s2 >= 11 & s2-s1 >=2) return('2')
  
  return (F)
}


play_game_hiho = function(p, h = T){
  s1 = 0
  s2 = 0
  
  while(T){
    
   
    
    out <- HiHo_win(s1,s2)
    
    if (out == '1') { return ('1')}
    if (out == '2') { return ('2')}
  
    if (play_point(p)){
      if ( h == T ){s1 = s1 + 1}
      else {h = T}
    }
    else {
      if ( h == F ) {s2 = s2 + 1}
      else {h = F}
    }
  
  }
  
}


play_game_par = function(p, h = T){
  s1 = 0
  s2 = 0
  
  while(T){
    
   
    
    out <- PAR_win(s1,s2)
    
    if (out == '1') {return ('1')}
    if (out == '2') {return ('2')}
    
    if (play_point(p)){s1 = s1 + 1}
    else {s2 = s2 + 1}
    }
    
}
  



play_match = function(p, h = T, game_fn){
  g1 = 0 
  g2 = 0
  
  while (T){
    if (g1 == 3) return('1')
    if (g2 == 3) return('2')
    
    if (game_fn(p, h) == '1') {g1 = g1 + 1; h = T}
    else {g2 = g2 + 1; h = F}
  }
}


res_hihoT = c()
res_hihoF = c()
res_par = c()
n = 1e4
for (p in seq(0, 0.5, 0.02)){
  t1 = c(); t2 = c(); t3 = c()
  
  for (i in 1:n) t1 = c(t1, play_game_hiho(p, h = T ))
  for (i in 1:n) t2 = c(t2, play_game_hiho(p, h = F ))
  for (i in 1:n) t3 = c(t3, play_game_par(p, h = F))
  
  res_hihoT = c(res_hihoT, sum(t1 =='1')/n)
  res_hihoF = c(res_hihoF, sum(t2 =='1')/n)
  res_par = c(res_par, sum(t3 =='1')/n)
  
}


par(mfrow=c(1,1))
plot(seq(0,0.5,0.02), res_hihoT, type = 'l', ylab = 'Probability of winning game', xlab = 'Probability of winning a single point')
lines(seq(0,0.5, 0.02), res_hihoF, col = 'red')
lines(seq(0,0.5, 0.02), res_par, col = 'blue')
abline(h=.5, lty = 2)
legend('topright', legend = c('English, serving', 'English, not serving', 'PAR'), cols = c('black', 'red', 'blue'), lty = 1)






f = function(v) c(v, 1-rev(v[-1]))






res_hihoT = c()
res_hihoF = c()
res_par = c()
n = 1e4
for (p in seq(0, 0.5, 0.02)){
  t1 = c(); t2 = c(); t3 = c()
  
  for (i in 1:n) t1 = c(t1, play_match(p, h = T, play_game_hiho ))
  for (i in 1:n) t2 = c(t2, play_match(p, h = F, play_game_hiho ))
  for (i in 1:n) t3 = c(t3, play_match(p, h = F, play_game_par ))
  
  res_hihoT = c(res_hihoT, sum(t1 =='1')/n)
  res_hihoF = c(res_hihoF, sum(t2 =='1')/n)
  res_par = c(res_par, sum(t3 =='1')/n)
  
}


par(mfrow=c(1,1))
plot(seq(0,1,0.02), f(res_hihoT), type = 'l', ylab = 'Probability of winning match', xlab = 'Probability of winning a single point')
lines(seq(0,1, 0.02), f(res_hihoF), col = 'red')
lines(seq(0,1, 0.02), f(res_par), col = 'blue')
abline(h=.5, lty = 2)
legend('topleft', legend = c('English, serving', 'English, not serving', 'PAR'), col = c('black', 'red', 'blue'), lty = 1)








for (i in 1:1e5) res = c(res, play_game(0.5, h = F, win_f = HiHo_win))

sum(res == '2')
