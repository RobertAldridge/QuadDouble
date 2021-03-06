
! tquaderq.f

module quadglobal
use qdmodule
implicit none
integer ndebug, ndigits, nerror, nquadl
end module

subroutine f_main

use qdmodule
use quadglobal
implicit none
integer i, kdebug, ndp, neps, nq1, nq2, n1
parameter (kdebug = 2, ndp = 64, neps = -64, nq1 = 9, nq2 = 8 * 2 ** nq1)
double precision dplog10q, d1, d2, second, tm0, tm1
type (qd_real) err, quaderq, fun01, fun02, fun03, fun04, fun05, fun06, fun07, &
  fun08, fun09, fun10, fun11, fun12, fun13, fun14, fun15a, fun15b, &
  t1, t2, t3, t4, wk(-1:nq2), xk(-1:nq2), x1, x2, gammax
external quaderq, fun01, fun02, fun03, fun04, fun05, fun06, fun07, fun08, &
  fun09, fun10, fun11, fun12, fun13, fun14, fun15a, fun15b, second, gammax

! integer*4 old_cw
! call f_fpu_fix_start (old_cw)

ndebug = kdebug
ndigits = ndp
nerror = 0
nquadl = nq1
write (6, 1) ndigits, neps, nquadl
1 format ('Quaderq test'/'Digits =',i6,'  Epsilon =',i6,'   Quadlevel =',i6)

tm0 = second ()
call initqerq (nq1, nq2, wk, xk)
tm1 = second ()
if(nerror > 0) stop
write (6, 2) tm1 - tm0
2 format ('Quadrature initialization completed: cpu time =',f12.6)

write (6, 11)
11 format (/'Continuous functions on finite itervals:'//&
  'Problem 1: Int_0^1 t*log(1+t) dt = 1/4')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quaderq (fun01, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
3 format ('Quadrature completed: CPU time =',f12.6/'Result =')
call qdwrite (6, t1)
t2 = 0.25d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1
4 format ('Actual error =',f10.6,'x10^',i5)

write (6, 12)
12 format (/'Problem 2: Int_0^1 t^2*arctan(t) dt = (pi - 2 + 2*log(2) )/12')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quaderq (fun02, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = (qdpi() - 2.d0 + 2.d0 * log (qdreal (2.d0) ) ) / 12.d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 13)
13 format (/'Problem 3: Int_0^(pi/2) e^t*cos(t) dt = 1/2*(e^(pi/2) - 1)')
x1 = 0.d0
x2 = 0.5d0 * qdpi()
tm0 = second ()
t1 = quaderq (fun03, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = 0.5d0 * (exp (0.5d0 * qdpi() ) - 1.d0)
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 14)
14 format (/ &
  'Problem 4: Int_0^1 arctan(sqrt(2+t^2) )/( (1+t^2)sqrt(2+t^2) ) dt = 5*Pi^2/96')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quaderq (fun04, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = 5.d0 * qdpi()**2 / 96.d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 15)
15 format (/&
  'Continuous functions on finite itervals, but non-diff at an endpoint'// &
  'Problem 5: Int_0^1 sqrt(t)*log(t) dt = -4/9')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quaderq (fun05, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = qdreal (-4.d0) / 9.d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 16)
16 format (/'Problem 6: Int_0^1 sqrt(1-t^2) dt = pi/4')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quaderq (fun06, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = 0.25d0 * qdpi()
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 17)
17 format (/&
  'Functions on finite intervals with integrable singularity at an endpoint.'//&
  'Problem 7: Int_0^1 sqrt(t)/sqrt(1-t^2) dt = 2*sqrt(pi)*gamma(3/4)/gamma(1/4)')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quaderq (fun07, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = 2.d0 * sqrt (qdpi() ) * gammax (qdreal (0.75d0) ) / gammax (qdreal (0.25d0) )
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 18)
18 format (/'Problem 8: Int_0^1 log(t)^2 dt = 2')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quaderq (fun08, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = 2.d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 19)
19 format (/'Problem 9: Int_0^(pi/2) log(cos(t) ) dt = -pi*log(2)/2')
x1 = 0.d0
x2 = 0.5d0 * qdpi()
tm0 = second ()
t1 = quaderq (fun09, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = -0.5d0 * qdpi() * log (qdreal (2.d0) )
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 20)
20 format (/'Problem 10: Int_0^(pi/2) sqrt(tan(t) ) dt = pi*sqrt(2)/2')
x1 = 0.d0
x2 = 0.5d0 * qdpi()
tm0 = second ()
t1 = quaderq (fun10, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = 0.5d0 * qdpi() * sqrt (qdreal (2.d0) )
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 21)
21 format (/&
  'Functions on an infinite interval (requiring a two-step solution'//&
  'Problem 11: Int_0^inf 1/(1+t^2) dt = pi/2')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quaderq (fun11, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = 0.5d0 * qdpi()
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 22)
22 format (/'Problem 12: Int_0^inf e^(-t)/sqrt(t) dt = sqrt(pi)')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quaderq (fun12, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = sqrt (qdpi() )
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 23)
23 format (/'Problem 13: Int_0^inf e^(-t^2/2) dt = sqrt(pi/2)')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quaderq (fun13, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = sqrt (0.5d0 * qdpi() )
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 24)
24 format (/&
  'Oscillatory functions on an infinite interval.'//&
  'Problem 14: Int_0^inf e^(-t)*cos(t) dt = 1/2')
x1 = 0.d0
x2 = 1.d0
tm0 = second ()
t1 = quaderq (fun14, x1, x2, nq1, nq2, wk, xk)
tm1 = second ()
write (6, 3) tm1 - tm0
call qdwrite (6, t1)
t2 = 0.5d0
call decmdq (t2 - t1, d1, n1)
write (6, 4) d1, n1

write (6, 25)
25 format (/'Problem 15: Int_0^inf sin(t)/t = pi/2')
x1 = 0.d0
x2 = qdpi()
tm0 = second ()
t1 = quaderq (fun15a, x1, x2, nq1, nq2, wk, xk)
x2 = 1.d0 / qdpi()
t2 = quaderq (fun15b, x1, x2, nq1, nq2, wk, xk)
t3 = t1 + 40320.d0 * t2 - 1.d0 / qdpi() + 2.d0 / qdpi() ** 3 &
  - 24.d0 / qdpi() ** 5 + 720.d0 / qdpi() ** 7
tm1 = second ()
write (6, 3) tm1 - tm0
t4 = 0.5d0 * qdpi()
call decmdq (t4 - t3, d1, n1)
write (6, 4) d1, n1
write (6, 26)
26 format ('Prob 15 error may be 40,000 X higher than estimated error.')

! call f_fpu_fix_end (old_cw)

stop
end

function fun01 (t)

use qdmodule
implicit none
type (qd_real) fun01, t

fun01 = t * log (1.d0 + t)
return
end

function fun02 (t)

use qdmodule
implicit none
type (qd_real) fun02, t

fun02 = t ** 2 * atan (t)
return
end

function fun03 (t)

use qdmodule
implicit none
type (qd_real) fun03, t

fun03 = exp(t) * cos(t)
return
end

function fun04 (t)

use qdmodule
implicit none
type (qd_real) fun04, t, t1

t1 = sqrt (2.d0 + t**2)
fun04 = atan(t1)/( (1.d0 + t**2)*t1)
return
end

function fun05 (t)

use qdmodule
implicit none
type (qd_real) fun05, t

fun05 = sqrt (t) * log (t)
return
end

function fun06 (t)

use qdmodule
implicit none
type (qd_real) fun06, t

fun06 = sqrt (1.d0 - t**2)
return
end

function fun07 (t)

use qdmodule
implicit none
type (qd_real) fun07, t

fun07 = sqrt(t) / sqrt (1.d0 - t**2)
return
end

function fun08 (t)

use qdmodule
implicit none
type (qd_real) fun08, t

fun08 = log (t) ** 2
return
end

function fun09 (t)

use qdmodule
implicit none
type (qd_real) fun09, t

fun09 = log (cos (t) )
return
end

function fun10 (t)

use qdmodule
implicit none
type (qd_real) fun10, t

fun10 = sqrt (tan (t) )
return
end

function fun11 (t)

use qdmodule
implicit none
type (qd_real) fun11, t

fun11 = 1.d0 / (1.d0 - 2.d0 * t + 2.d0 * t ** 2)
return
end

function fun12 (t)

use qdmodule
implicit none
type (qd_real) fun12, t, t1

t1 = 1.d0 / t - 1.d0
fun12 = exp (-t1) / sqrt (t ** 3 - t ** 4)
return
end

function fun13 (t)

use qdmodule
implicit none
type (qd_real) fun13, t, t1

t1 = 1.d0 / t - 1.d0
fun13 = exp (-0.5d0 * t1 ** 2) / t ** 2
return
end

function fun14 (t)

use qdmodule
implicit none
type (qd_real) fun14, t, t1

t1 = 1.d0 / t - 1.d0
fun14 = exp (-t1) * cos (t1) / t ** 2
return
end

function fun15a (t)

use qdmodule
use quadglobal
implicit none
type (qd_real) fun15a, t

fun15a = sin (t) / t
return
end

function fun15b (t)

use qdmodule
use quadglobal
implicit none
type (qd_real) fun15b, t

if(abs (t) > 1.d-10) then
  fun15b = t**7 * sin (1.d0 / t)
else
  fun15b = 0.d0
endif
return
end

subroutine initqerq (nq1, nq2, wk, xk)

use qdmodule
use quadglobal
implicit none
integer i, ierror, iprint, j, k, k1, nq1, nq2, ntab, ntabx
real*8 eps, h
parameter (iprint = 1000, ntabx = 1000)
type (qd_real) erfc, etab(ntabx), p2, spi, t1, t2, t3, t4, t5, &
  wk(-1:nq2), xk(-1:nq2)
external erfc

if(ndebug >= 1) then
  write (6, 1)
1 format ('initqerq: Error function quadrature initialization')
endif

eps = 1.d-64
p2 = 0.5d0 * qdpi()
spi = 2.d0 / sqrt (qdpi() )
h = 0.5d0 ** (nq1 - 2)
wk(-1) = dble (nq1)

do k = 0, nq2
  if(ndebug >= 2 .and. mod (k, iprint) == 0) write (6, *) k, nq2
  t1 = dble (k) * h
  xk(k) = erfc (t1, ntab, ntabx, etab)
  wk(k) = spi * exp (- t1 ** 2)
  if(wk(k) < eps) goto 100
enddo

write (6, 2) nq2
2 format ('initqerq: Table space parameter is too small; value =',i8)
nerror = 91
goto 130

100 continue

xk(-1) = dble (k)
if(ndebug >= 2) then
  write (6, 3) k
3 format ('initqerq: Table spaced used =',i8)
endif
goto 130

120 continue

nerror = ierror + 100
write (6, 4) nerror
4 format ('initqerq: Error in quadrature initialization; code =',i5)

130  continue

return
end

function quaderq (fun, x1, x2, nq1, nq2, wk, xk)

use qdmodule
use quadglobal
implicit none
integer i, ierror, ip(0:100), iz1, iz2, izx, j, k, k1, k2, n, nds, nq1, nq2, &
  nqq1
parameter (izx = 4)
logical log1, log2
real*8 d1, d2, d3, d4, dplog10q, h
type (qd_real) ax, bx, c10, quaderq, eps, eps1, eps2, err, fun, &
  tsum, s1, s2, s3, t1, t2, t3, t4, tw1, tw2, twi1, twi2, twmx, &
  wk(-1:nq2), xk(-1:nq2), x1, x2, xki, xt1, xx1, xx2
external fun, dplog10q

ax = 0.5d0 * (x2 - x1)
bx = 0.5d0 * (x2 + x1)
tsum = 0.d0
s1 = 0.d0
s2 = 0.d0
h = 4.d0
c10 = 10.d0
eps = 1.d-64

if(wk(-1) < dble (nq1) ) then
  write (6, 1) nq1
1 format ('quaderq: quadrature arrays have not been initialized; nq1 =',i6)
  nerror = 70
  goto 140
endif
nqq1 = dble (wk(-1) )
n = dble (xk(-1) )

do k = 0, nqq1
  ip(k) = 2 ** k
enddo

do k = 1, nq1
  h = 0.5d0 * h
  s3 = s2
  s2 = s1
  k1 = ip(nqq1-k)
  k2 = ip(nqq1-k+1)
  iz1 = 0
  iz2 = 0
  twmx = 0.d0

  do i = 0, n, k1
    if(mod (i, k2) /= 0 .or. k == 1) then
      xki = xk(i)
      xt1 = 1.d0 - xki
      xx1 = - ax * xt1 + bx
      xx2 = ax * xt1 + bx
      log1 = xx1 > x1
      log2 = xx2 < x2

      if(log1 .and. iz1 < izx) then
        t1 = fun (xx1)
        tw1 = t1 * wk(i)
        twi1 = abs (tw1)
        if(twi1 < eps) then
          iz1 = iz1 + 1
        else
          iz1 = 0
        endif
      else
        t1 = 0.d0
        tw1 = 0.d0
      endif

      if(i > 0 .and. log2 .and. iz2 < izx) then
        t2 = fun (xx2)
        tw2 = t2 * wk(i)
        twi2 = abs (tw2)
        if(twi2 < eps) then
          iz2 = iz2 + 1
        else
          iz2 = 0
        endif
      else
        t2 = 0.d0
        tw2 = 0.d0
      endif

      tsum = tsum + tw1 + tw2
      twmx = max (twmx, abs (tw1), abs (tw2) )
    endif
  enddo

  s1 =  ax * h * tsum
  eps1 = twmx * eps
  eps2 = max (twi1, twi2)
  d1 = dplog10q (abs (s1 - s2) )
  d2 = dplog10q (abs (s1 - s3) )
  d3 = dplog10q (eps1) - 1.d0
  d4 = dplog10q (eps2) - 1.d0

  if(k <= 2) then
    err = 1.d0
  elseif(d1 .eq. -9999.d0) then
    err = 0.d0
  else
    err = c10 ** nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3, d4) ) )
  endif

  if(ndebug >= 2) then
    write (6, 2) k, nq1, nint (dplog10q (abs (err) ) )
2   format ('quaderq: Iteration',i3,' of',i3,'; est error = 10^',i5, &
      '; approx value =')
    call qdwrite (6, s1)
  endif
  if(k >= 3 .and. err < eps1) goto 140
  if(k >= 3 .and. err < eps2) goto 120
enddo

write (6, 3) nint (dplog10q (abs (err) ) ), nquadl
3 format ('quaderq: Estimated error = 10^',i5/&
  'Increase Quadlevel for greater accuracy. Current Quadlevel =',i4)
goto 140

120 continue

write (6, 4) nint (dplog10q (abs (err) ) ), ndigits
4 format ('quaderq: Estimated error = 10^',i5/&
  'Increase working prec (Digits) for greater accuracy. Current Digits =',i4)
goto 140

130 continue

if(ierror > 0) nerror = ierror + 100
write (6, 5) nerror
5 format ('quaderq: Error in quadrature calculation; code =',i5)
s1 = 0.d0

140 continue

quaderq = s1
return
end

function erfc (t, ntab, ntabx, etab)

use qdmodule
implicit none
integer i, j, k, n, ndp, ntab, ntabx
real*8 alpha, d1, d2, dpi, dlog10, dlog2, dplog10q, eps
type (qd_real) erfc, etab(ntabx), t, t1, t2, t3, t4, t5
type (qd_real) t6, t7, t8
save alpha, dlog10

ndp = 64
eps = 1.d-64
if(ntab == 0) then

  dpi = acos (-1.d0)
  dlog10 = log (10.d0)
  dlog2 = log (2.d0)
  d1 = dpi / sqrt (ndp * dlog10)
  ntab = ndp * dlog10 / dpi
  if(ntab > ntabx) then
    write (6, *) 'ntabx must be at least', ntab
    stop
  endif
  n = abs (int (log (d1) / dlog2) ) + 1
  alpha = 0.5d0 ** (n + 6) * anint (d1 * 2.d0 ** (n + 6) )

  t1 = - alpha ** 2
  t2 = exp (t1)
  t3 = t2 ** 2
  t4 = 1.d0

  do i = 1, ntab
    t4 = t2 * t4
    etab(i) = t4
    t2 = t2 * t3
  enddo
endif

if(t == 0.d0) then
  erfc = 1.d0
  goto 200
endif
t1 = 0.d0
t2 = t ** 2
t3 = exp (-t2)
t4 = eps / t3 * 1.d-4

do k = 1, ntab
  t5 = etab(k) / (k ** 2 * alpha ** 2 + t2)
  t1 = t1 + t5
  if(abs (t5) < t4) goto 110
enddo

110 continue

erfc = t3 * alpha * t / qdpi() * (1.d0 / t2 + 2.d0 * t1) &
       + 2.d0 / (1.d0 - exp (2.d0 * qdpi() * t / alpha) )

200 continue
return
end

function gammax (t)

use qdmodule
implicit none
integer i, j, k, ndp, neps, nt, nwords
double precision alpha, con1, con2, d1, d2
parameter (con1 = 1.151292547d0, con2 = 1.974476770d0)
type (qd_real) eps, gammax, sum1, sum2, t, t1, t2, t3, t4, z

neps = -64
ndp = 64
eps = 1d-64

if(abs (t) > 170.d0) then
  write (6, *) 'gamma: argument too large'
  goto 120
elseif(t == anint (t) ) then
  if(t <= 0.d0) then
    write (6, *) 'gamma: invalid negative argument'
    z = 0.d0
    goto 120
  endif
  nt = dble (t)
  t1 = 1.d0

  do i = 2, nt - 1
    t1 = dble (i) * t1
  enddo

  z = t1
  goto 120
endif

alpha = aint (con1 * ndp + 1.d0)
t1 = t
d2 = 0.25d0 * alpha**2
t3 = 1.d0 / t1
sum1 = t3

do j = 1, 1000000000
  t3 = t3 * d2 / (dble (j) * (t1 + dble (j) ) )
  sum1 = sum1 + t3
  if(abs (t3) < abs (sum1) * eps) goto 100
enddo

write (6, *) 'gamma: loop overflow 1'
sum1 = 0.d0

100 continue

sum1 = t1 * (0.5d0 * alpha) ** t1 * sum1
t1 = -t
t3 = 1.d0 / t1
sum2 = t3

do j = 1, 1000000000
  t3 = t3 * d2 / (dble (j) * (t1 + dble (j) ) )
  sum2 = sum2 + t3
  if(abs (t3) < abs (sum2) * eps) goto 110
enddo

write (6, *) 'gamma: loop overflow 2'
sum2 = 0.d0

110 continue

sum2 = t1 * (0.5d0 * alpha) ** t1 * sum2

z = sqrt (qdpi() * sum1 / (t * sin (qdpi() * t) * sum2) )

120 continue

gammax = z
return
end

function dplog10q (a)

use qdmodule
implicit none
integer ia
double precision da, dplog10q, t1
type (qd_real) a

da = a
if(da .eq. 0.d0) then
  dplog10q = -9999.d0
else
  dplog10q = log10 (abs (da) )
endif

100 continue
return
end

subroutine decmdq (a, b, ib)

use qdmodule
implicit none
integer ia, ib
double precision da, b, t1, xlt
parameter (xlt = 0.3010299956639812d0)
type (qd_real) a

da = a
if(da .ne. 0.d0) then
  t1 = log10 (abs (da) )
  ib = t1
  if(t1 .lt. 0.d0) ib = ib - 1
  b = sign (10.d0 ** (t1 - ib), da)
else
  b = 0.d0
  ib = 0
endif

return
end
