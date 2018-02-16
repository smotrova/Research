!
!     dfimslMD.F90 - Declare double-precision DF IMSL MATH/LIBRARY routines
!
!     Copyright (c) 1997 by Visual Numerics, Inc.  All Rights Reserved.
!

module dfimslmd	
!dec$objcomment lib:'imsl.lib'
!dec$objcomment lib:'imsls_err.lib'
!dec$objcomment lib:'imslmpistub.lib'
      use dfimslc	! IMSL common routines
      use dfimslmc	! IMSL common math routines
      use dfimslcd	! IMSL common double-precision routines

!
!     Chapter 1:  Linear Systems
!

      interface
        subroutine dlsarg (n, a, lda, b, ipath, x)
          integer    n, lda, ipath
          double precision a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2arg (n, a, lda, b, ipath, x, fac, ipvt,wk)
          integer    n, lda, ipath, ipvt(*)
          double precision a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlfcrg (n, a, lda, fac, ldfac, ipvt,rcond)
          integer    n, lda, ldfac, ipvt(*)
          double precision rcond, a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dl2crg (n, a, lda, fac, ldfac, ipvt,rcond, z)
          integer    n, lda, ldfac, ipvt(*)
          double precision rcond, a(lda,*), fac(ldfac,*), z(*)
        end subroutine
      end interface

      interface
        subroutine dlfirg (n, a, lda, fac, ldfac, ipvt, b,ipath, x, res)
          integer    n, lda, ldfac, ipath, ipvt(*)
          double precision a(lda,*), fac(ldfac,*), b(*), x(*), res(*)
        end subroutine
      end interface

      interface
        subroutine dlsacg (n, a, lda, b, ipath, x)
          integer    n, lda, ipath
          complex    *16 a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2acg (n, a, lda, b, ipath, x, fac, ipvt,wk)
          integer    n, lda, ipath, ipvt(*)
          complex    *16 a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlslcg (n, a, lda, b, ipath, x)
          integer    n, lda, ipath
          complex    *16 a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2lcg (n, a, lda, b, ipath, x, fac, ipvt,wk)
          integer    n, lda, ipath, ipvt(*)
          complex    *16 a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlfccg (n, a, lda, fac, ldfac, ipvt,rcond)
          integer    n, lda, ldfac, ipvt(*)
          double precision rcond
          complex    *16 a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dl2ccg (n, a, lda, fac, ldfac, ipvt,rcond, z)
          integer    n, lda, ldfac, ipvt(*)
          double precision rcond
          complex    *16 a(lda,*), fac(ldfac,*), z(*)
        end subroutine
      end interface

      interface
        subroutine dlftcg (n, a, lda, fac, ldfac, ipvt)
          integer    n, lda, ldfac, ipvt(*)
          complex    *16 a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dl2tcg (n, a, lda, fac, ldfac, ipvt,scale)
          integer    n, lda, ldfac, ipvt(*)
          complex    *16 a(lda,*), fac(ldfac,*), scale(*)
        end subroutine
      end interface

      interface
        subroutine dlfscg (n, fac, ldfac, ipvt, b, ipath, x)
          integer    n, ldfac, ipath, ipvt(*)
          complex    *16 fac(ldfac,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dlficg (n, a, lda, fac, ldfac, ipvt, b,ipath, x, res)
          integer    n, lda, ldfac, ipath, ipvt(*)
          complex    *16 a(lda,*), fac(ldfac,*), b(*), x(*), res(*)
        end subroutine
      end interface

      interface
        subroutine dlfdcg (n, fac, ldfac, ipvt, det1, det2)
          integer    n, ldfac, ipvt(*)
          double precision det2
          complex    *16 det1, fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dlincg (n, a, lda, ainv, ldainv)
          integer    n, lda, ldainv
          complex    *16 a(lda,*), ainv(ldainv,*)
        end subroutine
      end interface

      interface
        subroutine dl2ncg (n, a, lda, ainv, ldainv, wk, iwk)
          integer    n, lda, ldainv, iwk(*)
          complex    *16 a(lda,*), ainv(ldainv,*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlslrt (n, a, lda, b, ipath, x)
          integer    n, lda, ipath
          double precision a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dlfcrt (n, a, lda, ipath, rcond)
          integer    n, lda, ipath
          double precision rcond, a(lda,*)
        end subroutine
      end interface

      interface
        subroutine dl2crt (n, a, lda, ipath, rcond, z)
          integer    n, lda, ipath
          double precision rcond, a(lda,*), z(*)
        end subroutine
      end interface

      interface
        subroutine dlfdrt (n, a, lda, det1, det2)
          integer    n, lda
          double precision det1, det2, a(lda,*)
        end subroutine
      end interface

      interface
        subroutine dlinrt (n, a, lda, ipath, ainv, ldainv)
          integer    n, lda, ipath, ldainv
          double precision a(lda,*), ainv(ldainv,*)
        end subroutine
      end interface

      interface
        subroutine dlslct (n, a, lda, b, ipath, x)
          integer    n, lda, ipath
          complex    *16 a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dlfcct (n, a, lda, ipath, rcond)
          integer    n, lda, ipath
          double precision rcond
          complex    *16 a(lda,*)
        end subroutine
      end interface

      interface
        subroutine dl2cct (n, a, lda, ipath, rcond, z)
          integer    n, lda, ipath
          double precision rcond
          complex    *16 a(lda,*), z(*)
        end subroutine
      end interface

      interface
        subroutine dlfdct (n, a, lda, det1, det2)
          integer    n, lda
          double precision det2
          complex    *16 det1, a(lda,*)
        end subroutine
      end interface

      interface
        subroutine dlinct (n, a, lda, ipath, ainv, ldainv)
          integer    n, lda, ipath, ldainv
          complex    *16 a(lda,*), ainv(ldainv,*)
        end subroutine
      end interface

      interface
        subroutine dlsads (n, a, lda, b, x)
          integer    n, lda
          double precision a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2ads (n, a, lda, b, x, fac, wk)
          integer    n, lda
          double precision a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlslds (n, a, lda, b, x)
          integer    n, lda
          double precision a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2lds (n, a, lda, b, x, fac, wk)
          integer    n, lda
          double precision a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlfcds (n, a, lda, fac, ldfac, rcond)
          integer    n, lda, ldfac
          double precision rcond, a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dl2cds (n, a, lda, fac, ldfac, rcond, z)
          integer    n, lda, ldfac
          double precision rcond, a(lda,*), fac(ldfac,*), z(*)
        end subroutine
      end interface

      interface
        subroutine dlftds (n, a, lda, fac, ldfac)
          integer    n, lda, ldfac
          double precision a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dlfsds (n, fac, ldfac, b, x)
          integer    n, ldfac
          double precision fac(ldfac,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dlfids (n, a, lda, fac, ldfac, b, x, res)
          integer    n, lda, ldfac
          double precision a(lda,*), fac(ldfac,*), b(*), x(*), res(*)
        end subroutine
      end interface

      interface
        subroutine dlfdds (n, fac, ldfac, det1, det2)
          integer    n, ldfac
          double precision det1, det2, fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dlinds (n, a, lda, ainv, ldainv)
          integer    n, lda, ldainv
          double precision a(lda,*), ainv(ldainv,*)
        end subroutine
      end interface

      interface
        subroutine dl2nds (n, a, lda, ainv, ldainv, wk)
          integer    n, lda, ldainv
          double precision a(lda,*), ainv(ldainv,*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlsasf (n, a, lda, b, x)
          integer    n, lda
          double precision a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2asf (n, a, lda, b, x, fac, ipvt, wk)
          integer    n, lda, ipvt(*)
          double precision a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlslsf (n, a, lda, b, x)
          integer    n, lda
          double precision a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2lsf (n, a, lda, b, x, fac, ipvt, wk)
          integer    n, lda, ipvt(*)
          double precision a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlfcsf (n, a, lda, fac, ldfac, ipvt,rcond)
          integer    n, lda, ldfac, ipvt(*)
          double precision rcond, a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dl2csf (n, a, lda, fac, ldfac, ipvt,rcond, z)
          integer    n, lda, ldfac, ipvt(*)
          double precision rcond, a(*), fac(ldfac,*), z(*)
        end subroutine
      end interface

      interface
        subroutine dlftsf (n, a, lda, fac, ldfac, ipvt)
          integer    n, lda, ldfac, ipvt(*)
          double precision a(*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dlfssf (n, fac, ldfac, ipvt, b, x)
          integer    n, ldfac, ipvt(*)
          double precision fac(ldfac,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dlfisf (n, a, lda, fac, ldfac, ipvt, b, x,res)
          integer    n, lda, ldfac, ipvt(*)
          double precision a(lda,*), fac(ldfac,*), b(*), x(*), res(*)
        end subroutine
      end interface

      interface
        subroutine dlfdsf (n, fac, ldfac, ipvt, det1, det2)
          integer    n, ldfac, ipvt(*)
          double precision det1, det2, fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dlsadh (n, a, lda, b, x)
          integer    n, lda
          complex    *16 a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2adh (n, a, lda, b, x, fac, wk)
          integer    n, lda
          complex    *16 a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlsldh (n, a, lda, b, x)
          integer    n, lda
          complex    *16 a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2ldh (n, a, lda, b, x, fac, wk)
          integer    n, lda
          complex    *16 a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlfcdh (n, a, lda, fac, ldfac, rcond)
          integer    n, lda, ldfac
          double precision rcond
          complex    *16 a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dl2cdh (n, a, lda, fac, ldfac, rcond, z)
          integer    n, lda, ldfac
          double precision rcond
          complex    *16 a(lda,*), fac(ldfac,*), z(*)
        end subroutine
      end interface

      interface
        subroutine dlftdh (n, a, lda, fac, ldfac)
          integer    n, lda, ldfac
          complex    *16 a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dlfsdh (n, fac, ldfac, b, x)
          integer    n, ldfac
          complex    *16 fac(ldfac,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dlfidh (n, a, lda, fac, ldfac, b, x, res)
          integer    n, lda, ldfac
          complex    *16 a(lda,*), fac(ldfac,*), b(*), x(*), res(*)
        end subroutine
      end interface

      interface
        subroutine dlfddh (n, fac, ldfac, det1, det2)
          integer    n, ldfac
          double precision det1, det2
          complex    *16 fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dlsahf (n, a, lda, b, x)
          integer    n, lda
          complex    *16 a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2ahf (n, a, lda, b, x, fac, ipvt, wk)
          integer    n, lda, ipvt(*)
          complex    *16 a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlslhf (n, a, lda, b, x)
          integer    n, lda
          complex    *16 a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2lhf (n, a, lda, b, x, fac, ipvt, wk)
          integer    n, lda, ipvt(*)
          complex    *16 a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlfchf (n, a, lda, fac, ldfac, ipvt,rcond)
          integer    n, lda, ldfac, ipvt(*)
          double precision rcond
          complex    *16 a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dl2chf (n, a, lda, fac, ldfac, ipvt,rcond, z)
          integer    n, lda, ldfac, ipvt(*)
          double precision rcond
          complex    *16 a(lda,*), fac(ldfac,*), z(*)
        end subroutine
      end interface

      interface
        subroutine dlfthf (n, a, lda, fac, ldfac, ipvt)
          integer    n, lda, ldfac, ipvt(*)
          complex    *16 a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dlfshf (n, fac, ldfac, ipvt, b, x)
          integer    n, ldfac, ipvt(*)
          complex    *16 fac(ldfac,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dlfihf (n, a, lda, fac, ldfac, ipvt, b, x,res)
          integer    n, lda, ldfac, ipvt(*)
          complex    *16 a(lda,*), fac(ldfac,*), b(*), x(*), res(*)
        end subroutine
      end interface

      interface
        subroutine dlfdhf (n, fac, ldfac, ipvt, det1, det2)
          integer    n, ldfac, ipvt(*)
          double precision det1, det2
          complex    *16 fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dlsltr (n, c, d, e, b)
          integer    n
          double precision c(*), d(*), e(*), b(*)
        end subroutine
      end interface

      interface
        subroutine dlslcr (n, c, a, b, ijob, y, u, ir, is)
          integer    n, ijob
          double precision c(*), a(*), b(*), y(*), u(*)
        end subroutine
      end interface

      interface
        subroutine dlsarb (n, a, lda, nlca, nuca, b, ipath,x)
          integer    n, lda, nlca, nuca, ipath
          double precision a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2arb (n, a, lda, nlca, nuca, b, ipath,x, fac, ipvt, wk)
          integer    n, lda, nlca, nuca, ipath, ipvt(*)
          double precision a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlslrb (n, a, lda, nlca, nuca, b, ipath,x)
          integer    n, lda, nlca, nuca, ipath
          double precision a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2lrb (n, a, lda, nlca, nuca, b, ipath,x, fac, ipvt, wk)
          integer    n, lda, nlca, nuca, ipath, ipvt(*)
          double precision a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlfcrb (n, a, lda, nlca, nuca, fac, ldfac,ipvt, rcond)
          integer    n, lda, nlca, nuca, ldfac, ipvt(*)
          double precision rcond, a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dl2crb (n, a, lda, nlca, nuca, fac, ldfac,ipvt, rcond, z)
          integer    n, lda, nlca, nuca, ldfac, ipvt(*)
          double precision rcond, a(lda,*), fac(ldfac,*), z(*)
        end subroutine
      end interface

      interface
        subroutine dlftrb (n, a, lda, nlca, nuca, fac, ldfac,ipvt)
          integer    n, lda, nlca, nuca, ldfac, ipvt(*)
          double precision a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dl2trb (n, a, lda, nlca, nuca, fac, ldfac,ipvt, scale)
          integer    n, lda, nlca, nuca, ldfac, ipvt(*)
          double precision a(lda,*), fac(ldfac,*), scale(*)
        end subroutine
      end interface

      interface
        subroutine dlfsrb (n, fac, ldfac, nlca, nuca, ipvt,b, ipath, x)
          integer    n, ldfac, nlca, nuca, ipath, ipvt(*)
          double precision fac(ldfac,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dlfirb (n, a, lda, nlca, nuca, fac, ldfac,ipvt, b, ipath, x, res)
          integer    n, lda, nlca, nuca, ldfac, ipath, ipvt(*)
          double precision a(lda,*), fac(ldfac,*), b(*), x(*), res(*)
        end subroutine
      end interface

      interface
        subroutine dlfdrb (n, fac, ldfac, nlca, nuca, ipvt,det1, det2)
          integer    n, ldfac, nlca, nuca, ipvt(*)
          double precision det1, det2, fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dlsaqs (n, a, lda, ncoda, b, x)
          integer    n, lda, ncoda
          double precision a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2aqs (n, a, lda, ncoda, b, x, fac, wk)
          integer    n, lda, ncoda
          double precision a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlslqs (n, a, lda, ncoda, b, x)
          integer    n, lda, ncoda
          double precision a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2lqs (n, a, lda, ncoda, b, x, fac, wk)
          integer    n, lda, ncoda
          double precision a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlslpb (n, a, lda, ncoda, ijob, u)
          integer    n, lda, ncoda, ijob
          double precision a(-ncoda+1:lda-ncoda,*), u(*)
        end subroutine
      end interface

      interface
        subroutine dl2lpb (n, a, lda, ncoda, ijob, u, w)
          integer    lda, n, ncoda, ijob
          double precision a(-ncoda+1:lda-ncoda,*), u(*), w(*)
        end subroutine
      end interface

      interface
        subroutine dlfcqs (n, a, lda, ncoda, fac, ldfac,rcond)
          integer    n, lda, ncoda, ldfac
          double precision rcond, a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dl2cqs (n, a, lda, ncoda, fac, ldfac,rcond, z)
          integer    n, lda, ncoda, ldfac
          double precision rcond, a(lda,*), fac(ldfac,*), z(*)
        end subroutine
      end interface

      interface
        subroutine dlftqs (n, a, lda, ncoda, fac, ldfac)
          integer    n, lda, ncoda, ldfac
          double precision a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dlfsqs (n, fac, ldfac, ncoda, b, x)
          integer    n, ldfac, ncoda
          double precision fac(ldfac,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dlfiqs (n, a, lda, ncoda, fac, ldfac, b,x, res)
          integer    n, lda, ncoda, ldfac
          double precision a(lda,*), fac(ldfac,*), b(*), x(*), res(*)
        end subroutine
      end interface

      interface
        subroutine dlfdqs (n, fac, ldfac, ncoda, det1, det2)
          integer    n, ldfac, ncoda
          double precision det1, det2, fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dlsltq (n, c, d, e, b)
          integer    n
          complex    *16 c(*), d(*), e(*), b(*)
        end subroutine
      end interface

      interface
        subroutine dlslcq (n, c, a, b, ijob, y, u, ir, is)
          integer    n, ijob
          double precision u(*)
          complex    *16 c(*), a(*), b(*), y(*)
        end subroutine
      end interface

      interface
        subroutine dlsacb (n, a, lda, nlca, nuca, b, ipath,x)
          integer    n, lda, nlca, nuca, ipath
          complex    *16 a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2acb (n, a, lda, nlca, nuca, b, ipath,x, fac, ipvt, wk)
          integer    n, lda, nlca, nuca, ipath, ipvt(*)
          complex    *16 a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlslcb (n, a, lda, nlca, nuca, b, ipath,x)
          integer    n, lda, nlca, nuca, ipath
          complex    *16 a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2lcb (n, a, lda, nlca, nuca, b, ipath,x, fac, ipvt, wk)
          integer    n, lda, nlca, nuca, ipath, ipvt(*)
          complex    *16 a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlfccb (n, a, lda, nlca, nuca, fac, ldfac,ipvt, rcond)
          integer    n, lda, nlca, nuca, ldfac, ipvt(*)
          double precision rcond
          complex    *16 a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dl2ccb (n, a, lda, nlca, nuca, fac, ldfac,ipvt, rcond, z)
          integer    n, lda, nlca, nuca, ldfac, ipvt(*)
          double precision rcond
          complex    *16 a(lda,*), fac(ldfac,*), z(*)
        end subroutine
      end interface

      interface
        subroutine dlftcb (n, a, lda, nlca, nuca, fac, ldfac,ipvt)
          integer    n, lda, nlca, nuca, ldfac, ipvt(*)
          complex    *16 a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dl2tcb (n, a, lda, nlca, nuca, fac, ldfac,ipvt, scale)
          integer    n, lda, nlca, nuca, ldfac, ipvt(*)
          complex    *16 a(lda,*), fac(ldfac,*), scale(*)
        end subroutine
      end interface

      interface
        subroutine dlfscb (n, fac, ldfac, nlca, nuca, ipvt,b, ipath, x)
          integer    n, ldfac, nlca, nuca, ipath, ipvt(*)
          complex    *16 fac(ldfac,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dlficb (n, a, lda, nlca, nuca, fac, ldfac,ipvt, b, ipath, x, res)
          integer    n, lda, nlca, nuca, ldfac, ipath, ipvt(*)
          complex    *16 a(lda,*), fac(ldfac,*), b(*), x(*), res(*)
        end subroutine
      end interface

      interface
        subroutine dlfdcb (n, fac, ldfac, nlca, nuca, ipvt,det1, det2)
          integer    n, ldfac, nlca, nuca, ipvt(*)
          double precision det2
          complex    *16 det1, fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dlsaqh (n, a, lda, ncoda, b, x)
          integer    n, lda, ncoda
          complex    *16 a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2aqh (n, a, lda, ncoda, b, x, fac, wk)
          integer    n, lda, ncoda
          complex    *16 a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlslqh (n, a, lda, ncoda, b, x)
          integer    n, lda, ncoda
          complex    *16 a(lda,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2lqh (n, a, lda, ncoda, b, x, fac, wk)
          integer    n, lda, ncoda
          complex    *16 a(lda,*), b(*), x(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlslqb (n, a, lda, ncoda, ijob, u)
          integer    n, lda, ncoda, ijob
          double precision a(-ncoda+1:lda-ncoda,*), u(*)
        end subroutine
      end interface

      interface
        subroutine dl2lqb (n, a, lda, ncoda, ijob, u, x, y)
          integer    lda, n, ncoda, ijob
          double precision a(-ncoda+1:lda-ncoda,*), u(*), x(*), y(*)
        end subroutine
      end interface

      interface
        subroutine dlfcqh (n, a, lda, ncoda, fac, ldfac,rcond)
          integer    n, lda, ncoda, ldfac
          double precision rcond
          complex    *16 a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dl2cqh (n, a, lda, ncoda, fac, ldfac,rcond, z)
          integer    n, lda, ncoda, ldfac
          double precision rcond
          complex    *16 a(lda,*), fac(ldfac,*), z(*)
        end subroutine
      end interface

      interface
        subroutine dlftqh (n, a, lda, ncoda, fac, ldfac)
          integer    n, lda, ncoda, ldfac
          complex    *16 a(lda,*), fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dlfsqh (n, fac, ldfac, ncoda, b, x)
          integer    n, ldfac, ncoda
          complex    *16 fac(ldfac,*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dlfiqh (n, a, lda, ncoda, fac, ldfac, b,x, res)
          integer    n, lda, ncoda, ldfac
          complex    *16 a(lda,*), fac(ldfac,*), b(*), x(*), res(*)
        end subroutine
      end interface

      interface
        subroutine dlfdqh (n, fac, ldfac, ncoda, det1, det2)
          integer    n, ldfac, ncoda
          double precision det1, det2
          complex    *16 fac(ldfac,*)
        end subroutine
      end interface

      interface
        subroutine dlslxg (n, nz, a, irow, jcol, b, ipath,iparam, rparam, x)
          integer    n, nz, ipath, irow(*), jcol(*), iparam(*)
          double precision a(*), b(*), rparam(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2lxg (n, nz, a, irow, jcol, b, ipath,iparam, rparam, x, wk, lwk, iwk, liwk)
          integer    n, nz, ipath, lwk, liwk, irow(*), jcol(*), iparam(*),  &
     &           iwk(*)
          double precision a(*), b(*), rparam(*), x(*), wk(lwk)
        end subroutine
      end interface

      interface
        subroutine dl4lxg (iparam, rparam)
          integer    iparam(*)
          double precision rparam(*)
        end subroutine
      end interface

      interface
        subroutine dlftxg (n, nz, a, irow, jcol, iparam,rparam, nfac, nl, fac, irfac, jcfac, ipvt, jpvt)
          integer    n, nz, nfac, nl, irow(*), jcol(*), iparam(*),          &
     &           irfac(*), jcfac(*), ipvt(*), jpvt(*)
          double precision a(*), rparam(*), fac(*)
        end subroutine
      end interface

      interface
        subroutine dl2txg (n, nz, a, irow, jcol, iparam,rparam, nfac, nl, fac, irfac, jcfac, ipvt, jpvt, wk,lwk, iwk, liwk)
          integer    n, nz, nfac, nl, lwk, liwk, irow(*), jcol(*),          &
     &           iparam(*), irfac(*), jcfac(*), ipvt(*), jpvt(*),       &
     &           iwk(*)
          double precision a(*), rparam(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlfsxg (n, nfac, nl, fac, irfac, jcfac,ipvt, jpvt, b, ipath, x)
          integer    n, nfac, nl, ipath, irfac(*), jcfac(*), ipvt(*),       &
     &           jpvt(*)
          double precision fac(*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dlslzg (n, nz, a, irow, jcol, b, ipath,iparam, rparam, x)
          integer    n, nz, ipath, irow(*), jcol(*), iparam(*)
          double precision rparam(*)
          complex    *16 a(*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2lzg (n, nz, a, irow, jcol, b, ipath,iparam, rparam, x, wk, lwk, iwk, liwk)
          integer    n, nz, ipath, lwk, liwk, irow(*), jcol(*), iparam(*),  &
     &           iwk(*)
          double precision rparam(*)
          complex    *16 a(*), b(*), x(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dl4lzg (iparam, rparam)
          integer    iparam(*)
          double precision rparam(*)
        end subroutine
      end interface

      interface
        subroutine dlftzg (n, nz, a, irow, jcol, iparam,rparam, nfac, nl, fac, irfac, jcfac, ipvt, jpvt)
          integer    n, nz, nfac, nl, irow(*), jcol(*), iparam(*),          &
     &           irfac(*), jcfac(*), ipvt(*), jpvt(*)
          double precision rparam(*)
          complex    *16 a(*), fac(*)
        end subroutine
      end interface

      interface
        subroutine dl2tzg (n, nz, a, irow, jcol, iparam,rparam, nfac, nl, fac, irfac, jcfac, ipvt, jpvt, wk,lwk, iwk, liwk)
          integer    n, nz, nfac, nl, lwk, liwk, irow(*), jcol(*),          &
     &           iparam(*), irfac(*), jcfac(*), ipvt(*), jpvt(*),       &
     &           iwk(*)
          double precision rparam(*)
          complex    *16 a(*), fac(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlfszg (n, nfac, nl, fac, irfac, jcfac,ipvt, jpvt, b, ipath, x)
          integer    n, nfac, nl, ipath, irfac(*), jcfac(*), ipvt(*),       &
     &           jpvt(*)
          complex    *16 fac(*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dlslxd (n, nz, a, irow, jcol, b, itwksp,x)
          integer    n, nz, itwksp, irow(*), jcol(*)
          double precision a(*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2lxd (n, nz, a, irow, jcol, b, x, iper,iparam, rparam, wk, lwk, iwk, liwk)
          integer    n, nz, lwk, liwk, irow(*), jcol(*), iper(*),           &
     &           iparam(*), iwk(*)
          double precision a(*), b(*), x(*), rparam(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dl4lxd (iparam, rparam)
          integer    iparam(*)
          double precision rparam(*)
        end subroutine
      end interface

      interface
        subroutine dlscxd (n, nz, irow, jcol, ijob, itwksp,maxsub, nzsub, inzsub, maxnz, ilnz, iper, invper,ispace)
          integer    n, nz, ijob, itwksp, maxsub, maxnz, ispace, irow(*),   &
     &           jcol(*), nzsub(*), inzsub(*), ilnz(*), iper(*),        &
     &           invper(*)
        end subroutine
      end interface

      interface
        subroutine dl2cxd (n, nz, irow, jcol, ijob, maxsub,nzsub, inzsub, maxnz, ilnz, iper, invper, ispace, liwk,iwk)
          integer    n, nz, ijob, maxsub, maxnz, ispace, liwk, irow(*),     &
     &           jcol(*), nzsub(*), inzsub(*), ilnz(*), iper(*),        &
     &           invper(*), iwk(*)
        end subroutine
      end interface

      interface
        subroutine dlnfxd (n, nz, a, irow, jcol, ijob,maxsub, nzsub, inzsub, maxnz, ilnz, iper, invper,ispace, itwksp, diag, rlnz, rparam)
          integer    n, nz, ijob, maxsub, maxnz, ispace, itwksp, irow(*),   &
     &           jcol(*), nzsub(*), inzsub(*), ilnz(*), iper(*),        &
     &           invper(*)
          double precision a(*), diag(*), rlnz(*), rparam(*)
        end subroutine
      end interface

      interface
        subroutine dl2fxd (n, nz, a, irow, jcol, ijob,maxsub, nzsub, inzsub, maxnz, ilnz, iper, invper,ispace, diag, rlnz, rparam, wk, lwk, iwk, liwk)
          integer    n, nz, ijob, maxsub, maxnz, ispace, lwk, liwk,         &
     &           irow(*), jcol(*), nzsub(*), inzsub(*), ilnz(*),        &
     &           iper(*), invper(*), iwk(*)
          double precision a(*), diag(*), rlnz(*), rparam(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlfsxd (n, maxsub, nzsub, inzsub, maxnz,rlnz, ilnz, diag, iper, b, x)
          integer    n, maxsub, maxnz, nzsub(*), inzsub(*), ilnz(*),        &
     &           iper(*)
          double precision rlnz(*), diag(*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dlslzd (n, nz, a, irow, jcol, b, itwksp,x)
          integer    n, nz, itwksp, irow(*), jcol(*)
          complex    *16 a(*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2lzd (n, nz, a, irow, jcol, b, x, iper,iparam, rparam, wk, lwk, iwk, liwk)
          integer    n, nz, lwk, liwk, irow(*), jcol(*), iper(*),           &
     &           iparam(*), iwk(*)
          double precision rparam(*)
          complex    *16 a(*), b(*), x(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dl4lzd (iparam, rparam)
          integer    iparam(*)
          double precision rparam(*)
        end subroutine
      end interface

      interface
        subroutine dlnfzd (n, nz, a, irow, jcol, ijob,maxsub, nzsub, inzsub, maxnz, ilnz, iper, invper,ispace, itwksp, diag, rlnz, rparam)
          integer    n, nz, ijob, maxsub, maxnz, ispace, itwksp, irow(*),   &
     &           jcol(*), nzsub(*), inzsub(*), ilnz(*), iper(*),        &
     &           invper(*)
          double precision rparam(*)
          complex    *16 a(*), diag(*), rlnz(*)
        end subroutine
      end interface

      interface
        subroutine dl2fzd (n, nz, a, irow, jcol, ijob,maxsub, nzsub, inzsub, maxnz, ilnz, iper, invper,ispace, diag, rlnz, rparam, wk, lwk, iwk, liwk)
          integer    n, nz, ijob, maxsub, maxnz, ispace, lwk, liwk,         &
     &           irow(*), jcol(*), nzsub(*), inzsub(*), ilnz(*),        &
     &           iper(*), invper(*), iwk(*)
          double precision rparam(*)
          complex    *16 a(*), diag(*), rlnz(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlfszd (n, maxsub, nzsub, inzsub, maxnz,rlnz, ilnz, diag, iper, b, x)
          integer    n, maxsub, maxnz, nzsub(*), inzsub(*), ilnz(*),        &
     &           iper(*)
          complex    *16 rlnz(*), diag(*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dlslto (n, a, b, ipath, x)
          integer    n, ipath
          double precision a(*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2lto (n, a, b, ipath, x, wk)
          integer    n, ipath
          double precision a(*), b(*), x(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlsltc (n, a, b, ipath, x)
          integer    n, ipath
          complex    *16 a(*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2ltc (n, a, b, ipath, x, wk)
          integer    n, ipath
          complex    *16 a(*), b(*), x(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlslcc (n, a, b, ipath, x)
          integer    n, ipath
          complex    *16 a(*), b(*), x(*)
        end subroutine
      end interface

      interface
        subroutine dl2lcc (n, a, b, ipath, x, ahat, wk)
          integer    n, ipath
          double precision wk(*)
          complex    *16 a(*), b(*), x(*), ahat(*)
        end subroutine
      end interface

      interface
        subroutine dpcgrc (ido, n, x, p, r, z, relerr, itmax)
          integer    ido, n, itmax
          double precision relerr, x(*), p(*), r(*), z(*)
        end subroutine
      end interface

      interface
        subroutine dp2grc (ido, n, x, p, r, z, relerr, itmax,tri, work, iwork)
          integer    ido, n, itmax, iwork(*)
          double precision relerr, x(*), p(*), r(*), z(*), tri(2,*),        &
     &           work(*)
        end subroutine
      end interface

      interface
        subroutine djcgrc (ido, n, diag, x, p, r, z, relerr,itmax)
          integer    ido, n, itmax
          double precision relerr, diag(*), x(*), p(*), r(*), z(*)
        end subroutine
      end interface

      interface
        subroutine dj2grc (ido, n, diag, x, p, r, z, relerr,itmax, tri, work, iwork)
          integer    ido, n, itmax, iwork(*)
          double precision relerr, diag(*), x(*), p(*), r(*), z(*),         &
     &           tri(2,*), work(*)
        end subroutine
      end interface

      interface
        subroutine dlqrrv (nra, nca, numexc, a, lda, x, ldx)
          integer    nra, nca, numexc, lda, ldx
          double precision a(lda,*), x(ldx,*)
        end subroutine
      end interface

      interface
        subroutine dl2rrv (nra, nca, numexc, a, lda, x, ldx,fac, ldfac, work)
          integer    nra, nca, numexc, lda, ldx, ldfac
          double precision a(lda,*), x(ldx,*), fac(ldfac,*),                &
     &           work(nca+numexc+1,*)
        end subroutine
      end interface

      interface
        subroutine dlclsq (nra, nca, ncon, a, lda, b, c, ldc,bl, bu, irtype, xlb, xub, x, res)
          integer    nra, nca, ncon, lda, ldc, irtype(*)
          double precision a(lda,*), b(*), c(ldc,*), bl(*), bu(*), xlb(*),  &
     &           xub(*), x(*), res(*)
        end subroutine
      end interface

      interface
        subroutine dl2lsq (nra, nca, ncon, a, lda, b, c, ldc,bl, bu, irtype, xlb, xub, x, res, wk, iwk)
          integer    nra, nca, ncon, lda, ldc, irtype(*), iwk(*)
          double precision a(lda,*), b(*), c(ldc,*), bl(*), bu(*), xlb(*),  &
     &           xub(*), x(*), res(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlqrrr (nra, nca, a, lda, pivot, ipvt, qr,ldqr, qraux, conorm)
          integer    nra, nca, lda, ldqr, ipvt(*)
          double precision a(*), qr(*), qraux(*), conorm(*)
          logical    pivot
        end subroutine
      end interface

      interface
        subroutine dl2rrr (nra, nca, a, lda, pivot, ipvt, qr,ldqr, qraux, conorm, work)
          integer    nra, nca, lda, ldqr, ipvt(*)
          double precision a(lda,*), qr(ldqr,*), qraux(*), conorm(*),       &
     &           work(*)
          logical    pivot
        end subroutine
      end interface

      interface
        subroutine dlqerr (nrqr, ncqr, qr, ldqr, qraux, q,ldq)
          integer    nrqr, ncqr, ldqr, ldq
          double precision qr(ldqr,*), qraux(*), q(ldq,*)
        end subroutine
      end interface

      interface
        subroutine dl2err (nrqr, ncqr, qr, ldqr, qraux, q,ldq, work)
          integer    nrqr, ncqr, ldqr, ldq
          double precision qr(ldqr,*), qraux(*), q(ldq,*), work(*)
        end subroutine
      end interface

      interface
        subroutine dlqrsl (nra, kbasis, qr, ldqr, qraux, b,ipath, qb, qtb, x, res, ab)
          integer    nra, kbasis, ldqr, ipath
          double precision qr(ldqr,*), qraux(*), b(*), qb(*), qtb(*),       &
     &           x(*), res(*), ab(*)
        end subroutine
      end interface

      interface
        subroutine dlupqr (nrow, ncol, alpha, w, y, q, ldq,r, ldr, ipath, qnew, ldqnew, rnew, ldrnew)
          integer    nrow, ncol, ldq, ldr, ipath, ldqnew, ldrnew
          double precision alpha, w(*), y(*), q(*), r(*), qnew(*), rnew(*)
        end subroutine
      end interface

      interface
        subroutine dl2pqr (nrow, ncol, alpha, w, y, q, ldq,r, ldr, ipath, qnew, ldqnew, rnew, ldrnew, z, work)
          integer    nrow, ncol, ldq, ldr, ipath, ldqnew, ldrnew
          double precision alpha, w(*), y(*), q(ldq,*), r(ldr,*),           &
     &           qnew(ldqnew,*), rnew(ldrnew,*), z(*), work(*)
        end subroutine
      end interface

      interface
        subroutine dlchrg (n, a, lda, pivot, ipvt, fac,ldfac)
          integer    n, lda, ldfac, ipvt(*)
          double precision a(lda,*), fac(ldfac,*)
          logical    pivot
        end subroutine
      end interface

      interface
        subroutine dlupch (n, r, ldr, x, rnew, ldrnew, cs,sn)
          integer    n, ldr, ldrnew
          double precision r(ldr,*), x(*), rnew(ldrnew,*), cs(*), sn(*)
        end subroutine
      end interface

      interface
        subroutine dldnch (n, r, ldr, x, rnew, ldrnew, cs,sn)
          integer    n, ldr, ldrnew
          double precision r(ldr,*), x(*), rnew(ldrnew,*), cs(*), sn(*)
        end subroutine
      end interface

      interface
        subroutine dlsvcr (nra, nca, a, lda, ipath, tol,irank, s, u, ldu, v, ldv)
          integer    nra, nca, lda, ipath, irank, ldu, ldv
          double precision tol
          complex    *16 a(lda,*), s(*), u(ldu,*), v(ldv,*)
        end subroutine
      end interface

      interface
        subroutine dl2vcr (nra, nca, a, lda, ipath, tol,irank, s, u, ldu, v, ldv, wka, wk)
          integer    nra, nca, lda, ipath, irank, ldu, ldv
          double precision tol
          complex    *16 a(lda,*), s(*), u(*), v(*), wka(*), wk(*)
        end subroutine
      end interface

!
!     Chapter 2:  Eigensystem Analysis
!

      interface
        subroutine devlrg (n, a, lda, eval)
          integer    n, lda
          double precision a(lda,*)
          complex    *16 eval(*)
        end subroutine
      end interface

      interface
        subroutine de3lrg (n, a, lda, eval, acopy, wk, iwk)
          integer    n, lda, iwk(n,2)
          double precision a(lda,*), acopy(*), wk(n,4)
          complex    *16 eval(*)
        end subroutine
      end interface

      interface
        subroutine devcrg (n, a, lda, eval, evec, ldevec)
          integer    n, lda, ldevec
          double precision a(lda,*)
          complex    *16 eval(*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine de8crg (n, a, lda, eval, evec, ldevec,acopy, ecopy, wk, iwk)
          integer    n, lda, ldevec, iwk(n)
          double precision a(lda,*), acopy(*), ecopy(*), wk(n,6)
          complex    *16 eval(*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        double precision function depirg (n, neval, a, lda,eval, evec, ldevec)
          integer    n, neval, lda, ldevec
          double precision a(lda,*)
          complex    *16 eval(*), evec(ldevec,*)
        end function
      end interface

      interface
        double precision function de2irg (n, neval, a, lda,eval, evec, ldevec, cwk)
          integer    n, neval, lda, ldevec
          double precision a(lda,*)
          complex    *16 eval(*), evec(ldevec,*), cwk(*)
        end function
      end interface

      interface
        subroutine devlcg (n, a, lda, eval)
          integer    n, lda
          complex    *16 a(lda,*), eval(*)
        end subroutine
      end interface

      interface
        subroutine de3lcg (n, a, lda, eval, acopy, rwk, cwk,iwk)
          integer    n, lda, iwk(n)
          double precision rwk(*)
          complex    *16 a(lda,*), eval(*), acopy(*), cwk(n,*)
        end subroutine
      end interface

      interface
        subroutine devccg (n, a, lda, eval, evec, ldevec)
          integer    n, lda, ldevec
          complex    *16 a(lda,*), eval(*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine de6ccg (n, a, lda, eval, evec, ldevec,acopy, rwk, cwk, iwk)
          integer    n, lda, ldevec, iwk(n)
          double precision rwk(n)
          complex    *16 a(lda,*), eval(*), evec(ldevec,*), acopy(*),       &
     &           cwk(n,*)
        end subroutine
      end interface

      interface
        double precision function depicg (n, neval, a, lda,eval, evec, ldevec)
          integer    n, neval, lda, ldevec
          complex    *16 a(lda,*), eval(*), evec(ldevec,*)
        end function
      end interface

      interface
        double precision function de2icg (n, neval, a, lda,eval, evec, ldevec, wk)
          integer    n, neval, lda, ldevec
          complex    *16 a(lda,*), eval(*), evec(ldevec,*), wk(*)
        end function
      end interface

      interface
        subroutine devlsf (n, a, lda, eval)
          integer    n, lda
          double precision a(lda,*), eval(*)
        end subroutine
      end interface

      interface
        subroutine de4lsf (n, a, lda, eval, wk, iwk)
          integer    n, lda, iwk(n)
          double precision a(lda,*), eval(*), wk(n,2)
        end subroutine
      end interface

      interface
        subroutine devasf (n, neval, a, lda, small, eval)
          integer    n, neval, lda
          double precision a(lda,*), eval(*)
          logical    small
        end subroutine
      end interface

      interface
        subroutine de4asf (n, neval, a, lda, small, eval, wk,iwk)
          integer    n, neval, lda, iwk(n)
          double precision a(lda,*), eval(*), wk(n,4)
          logical    small
        end subroutine
      end interface

      interface
        subroutine devesf (n, nevec, a, lda, small, eval,evec, ldevec)
          integer    n, nevec, lda, ldevec
          double precision a(lda,*), eval(*), evec(ldevec,*)
          logical    small
        end subroutine
      end interface

      interface
        subroutine de5esf (n, nevec, a, lda, small, eval,evec, ldevec, wk, iwk)
          integer    n, nevec, lda, ldevec, iwk(n)
          double precision a(lda,*), eval(*), evec(ldevec,*), wk(n,9)
          logical    small
        end subroutine
      end interface

      interface
        subroutine devbsf (n, mxeval, a, lda, elow, ehigh,neval, eval)
          integer    n, mxeval, lda, neval
          double precision elow, ehigh, a(lda,*), eval(mxeval)
        end subroutine
      end interface

      interface
        subroutine de5bsf (n, mxeval, a, lda, elow, ehigh,neval, eval, wk, iwk)
          integer    n, mxeval, lda, neval, iwk(n)
          double precision elow, ehigh, a(lda,*), eval(mxeval), wk(n,5)
        end subroutine
      end interface

      interface
        subroutine devfsf (n, mxeval, a, lda, elow, ehigh,neval, eval, evec, ldevec)
          integer    n, mxeval, lda, neval, ldevec
          double precision elow, ehigh, a(lda,*), eval(*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine de3fsf (n, mxeval, a, lda, elow, ehigh,neval, eval, evec, ldevec, wk, iwk)
          integer    n, mxeval, lda, neval, ldevec, iwk(n)
          double precision elow, ehigh, a(lda,*), eval(mxeval),             &
     &           evec(ldevec,*), wk(n,9)
        end subroutine
      end interface

      interface
        subroutine devlsb (n, a, lda, ncoda, eval)
          integer    n, lda, ncoda
          double precision a(lda,*), eval(*)
        end subroutine
      end interface

      interface
        subroutine de3lsb (n, a, lda, ncoda, eval, acopy, wk)
          integer    n, lda, ncoda
          double precision a(lda,*), eval(*), acopy(ncoda+1,*), wk(*)
        end subroutine
      end interface

      interface
        subroutine devcsb (n, a, lda, ncoda, eval, evec,ldevec)
          integer    n, lda, ncoda, ldevec
          double precision a(lda,*), eval(*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine de4csb (n, a, lda, ncoda, eval, evec,ldevec, acopy, wk, iwk)
          integer    n, lda, ncoda, ldevec, iwk(*)
          double precision a(lda,*), eval(*), evec(ldevec,*),               &
     &           acopy(ncoda+1,*), wk(*)
        end subroutine
      end interface

      interface
        subroutine devasb (n, neval, a, lda, ncoda, small,eval)
          integer    n, neval, lda, ncoda
          double precision a(lda,*), eval(*)
          logical    small
        end subroutine
      end interface

      interface
        subroutine de3asb (n, neval, a, lda, ncoda, small,eval, acopy, wk)
          integer    n, neval, lda, ncoda
          double precision a(lda,*), eval(*), acopy(ncoda+1,*), wk(3*n)
          logical    small
        end subroutine
      end interface

      interface
        subroutine devesb (n, nevec, a, lda, ncoda, small,eval, evec, ldevec)
          integer    n, nevec, lda, ncoda, ldevec
          double precision a(lda,*), eval(*), evec(ldevec,*)
          logical    small
        end subroutine
      end interface

      interface
        subroutine de4esb (n, nevec, a, lda, ncoda, small,eval, evec, ldevec, acopy, wk, iwk)
          integer    n, nevec, lda, ncoda, ldevec, iwk(n)
          double precision a(lda,*), eval(*), evec(ldevec,*),               &
     &           acopy(ncoda+1,*), wk(n*(2*ncoda+5))
          logical    small
        end subroutine
      end interface

      interface
        subroutine devbsb (n, mxeval, a, lda, ncoda, elow,ehigh, neval, eval)
          integer    n, mxeval, lda, ncoda, neval
          double precision elow, ehigh, a(lda,*), eval(*)
        end subroutine
      end interface

      interface
        subroutine de3bsb (n, mxeval, a, lda, ncoda, elow,ehigh, neval, eval, acopy, wk)
          integer    n, mxeval, lda, ncoda, neval
          double precision elow, ehigh, a(lda,*), eval(*),                  &
     &           acopy(ncoda+1,*), wk(n,*)
        end subroutine
      end interface

      interface
        subroutine devfsb (n, mxeval, a, lda, ncoda, elow,ehigh, neval, eval, evec, ldevec)
          integer    n, mxeval, lda, ncoda, neval, ldevec
          double precision elow, ehigh, a(lda,*), eval(*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine de3fsb (n, mxeval, a, lda, ncoda, elow,ehigh, neval, eval, evec, ldevec, acopy, wk1, wk2, iwk)
          integer    n, mxeval, lda, ncoda, neval, ldevec, iwk(n)
          double precision elow, ehigh, a(lda,*), eval(*), evec(ldevec,*),  &
     &           acopy(ncoda+1,*), wk1(6*n), wk2(2*n*ncoda+n)
        end subroutine
      end interface

      interface
        double precision function depisb (n, neval, a, lda,ncoda, eval, evec, ldevec)
          integer    n, neval, lda, ncoda, ldevec
          double precision a(lda,*), eval(*), evec(ldevec,*)
        end function
      end interface

      interface
        double precision function de2isb (n, neval, a, lda,ncoda, eval, evec, ldevec, wk)
          integer    n, neval, lda, ncoda, ldevec
          double precision a(lda,*), eval(*), evec(ldevec,*), wk(*)
        end function
      end interface

      interface
        subroutine devlhf (n, a, lda, eval)
          integer    n, lda
          double precision eval(*)
          complex    *16 a(lda,*)
        end subroutine
      end interface

      interface
        subroutine de3lhf (n, a, lda, eval, acopy, rwk, cwk,iwk)
          integer    n, lda, iwk(n)
          double precision eval(*), rwk(*)
          complex    *16 a(lda,*), acopy(*), cwk(*)
        end subroutine
      end interface

      interface
        subroutine devchf (n, a, lda, eval, evec, ldevec)
          integer    n, lda, ldevec
          double precision eval(*)
          complex    *16 a(lda,*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine de5chf (n, a, lda, eval, evec, ldevec,acopy, rwk, cwk, iwk)
          integer    n, lda, ldevec, iwk(n)
          double precision eval(*), rwk(*)
          complex    *16 a(lda,*), evec(ldevec,*), acopy(*), cwk(*)
        end subroutine
      end interface

      interface
        subroutine devahf (n, neval, a, lda, small, eval)
          integer    n, neval, lda
          double precision eval(*)
          complex    *16 a(lda,*)
          logical    small
        end subroutine
      end interface

      interface
        subroutine de3ahf (n, neval, a, lda, small, eval,acopy, rwk, cwk, iwk)
          integer    n, neval, lda, iwk(*)
          double precision eval(*), rwk(n,*)
          complex    *16 a(lda,*), acopy(n,*), cwk(n,*)
          logical    small
        end subroutine
      end interface

      interface
        subroutine devehf (n, nevec, a, lda, small, eval,evec, ldevec)
          integer    n, nevec, lda, ldevec
          double precision eval(*)
          complex    *16 a(lda,*), evec(ldevec,*)
          logical    small
        end subroutine
      end interface

      interface
        subroutine de3ehf (n, nevec, a, lda, small, eval,evec, ldevec, acopy, rw1, rw2, cwk, iwk)
          integer    n, nevec, lda, ldevec, iwk(*)
          double precision eval(*), rw1(n,*), rw2(n,*)
          complex    *16 a(lda,*), evec(ldevec,*), acopy(n,*), cwk(n,*)
          logical    small
        end subroutine
      end interface

      interface
        subroutine devbhf (n, mxeval, a, lda, elow, ehigh,neval, eval)
          integer    n, mxeval, lda, neval
          double precision elow, ehigh, eval(*)
          complex    *16 a(lda,*)
        end subroutine
      end interface

      interface
        subroutine de3bhf (n, mxeval, a, lda, elow, ehigh,neval, eval, acopy, rwk, cwk, iwk)
          integer    n, mxeval, lda, neval, iwk(*)
          double precision elow, ehigh, eval(*), rwk(n,*)
          complex    *16 a(lda,*), acopy(n,*), cwk(n,*)
        end subroutine
      end interface

      interface
        subroutine devfhf (n, mxeval, a, lda, elow, ehigh,neval, eval, evec, ldevec)
          integer    n, mxeval, lda, neval, ldevec
          double precision elow, ehigh, eval(*)
          complex    *16 a(lda,*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine de3fhf (n, mxeval, a, lda, elow, ehigh,neval, eval, evec, ldevec, acopy, ecopy, rwk, cwk, iwk)
          integer    n, mxeval, lda, neval, ldevec, iwk(*)
          double precision elow, ehigh, eval(*), ecopy(n,*), rwk(n,*)
          complex    *16 a(lda,*), evec(ldevec,*), acopy(n,*), cwk(n,*)
        end subroutine
      end interface

      interface
        double precision function depihf (n, neval, a, lda,eval, evec, ldevec)
          integer    n, neval, lda, ldevec
          double precision eval(*)
          complex    *16 a(lda,*), evec(ldevec,*)
        end function
      end interface

      interface
        double precision function de2ihf (n, neval, a, lda,eval, evec, ldevec, wk)
          integer    n, neval, lda, ldevec
          double precision eval(*)
          complex    *16 a(lda,*), evec(ldevec,*), wk(*)
        end function
      end interface

      interface
        subroutine devlrh (n, a, lda, eval)
          integer    n, lda
          double precision a(lda,*)
          complex    *16 eval(*)
        end subroutine
      end interface

      interface
        subroutine de3lrh (n, a, lda, eval, acopy, wk, iwk)
          integer    n, lda, iwk(n,1)
          double precision a(lda,*), acopy(n,*), wk(n,3)
          complex    *16 eval(*)
        end subroutine
      end interface

      interface
        subroutine devcrh (n, a, lda, eval, evec, ldevec)
          integer    n, lda, ldevec
          double precision a(lda,*)
          complex    *16 eval(*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine de6crh (n, a, lda, eval, evec, ldevec,acopy, ecopy, rwk, iwk)
          integer    n, lda, ldevec, iwk(n)
          double precision a(lda,*), acopy(n,*), ecopy(n,n), rwk(n,3)
          complex    *16 eval(*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine devlch (n, a, lda, eval)
          integer    n, lda
          complex    *16 a(lda,*), eval(*)
        end subroutine
      end interface

      interface
        subroutine de3lch (n, a, lda, eval, acopy, rwk, iwk)
          integer    n, lda, iwk(n)
          complex    *16 a(lda,*), eval(*), acopy(n,*)
          double precision rwk(n)
        end subroutine
      end interface

      interface
        subroutine devcch (n, a, lda, eval, evec, ldevec)
          integer    n, lda, ldevec
          complex    *16 a(lda,*), eval(*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine de4cch (n, a, lda, eval, evec, ldevec,acopy, cwork, rwk, iwk)
          integer    n, lda, ldevec, iwk(n)
          double precision rwk(n)
          complex    *16 a(lda,*), eval(*), evec(ldevec,*), acopy(n,*),     &
     &           cwork(n,*)
        end subroutine
      end interface

      interface
        subroutine dgvlrg (n, a, lda, b, ldb, alpha, beta)
          integer    n, lda, ldb
          double precision a(lda,*), b(ldb,*), beta(*)
          complex    *16 alpha(*)
        end subroutine
      end interface

      interface
        subroutine dg3lrg (n, a, lda, b, ldb, alpha, beta,acopy, bcopy, rwk, cwk, iwk)
          integer    n, lda, ldb, iwk(n)
          double precision a(lda,*), b(ldb,*), beta(*), acopy(*),           &
     &           bcopy(*), rwk(*)
          complex    *16 alpha(*), cwk(*)
        end subroutine
      end interface

      interface
        subroutine dgvcrg (n, a, lda, b, ldb, alpha, beta,evec, ldevec)
          integer    n, lda, ldb, ldevec
          double precision a(lda,*), b(ldb,*), beta(*)
          complex    *16 alpha(*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine dg8crg (n, a, lda, b, ldb, alpha, beta,evec, ldevec, acopy, bcopy, ecopy, rwk, cwk, iwk)
          integer    n, lda, ldb, ldevec, iwk(n)
          double precision a(lda,*), b(ldb,*), beta(*), acopy(*),           &
     &           bcopy(*), ecopy(*), rwk(*)
          complex    *16 alpha(*), evec(ldevec,*), cwk(*)
        end subroutine
      end interface

      interface
        double precision function dgpirg (n, neval, a, lda,b, ldb, alpha, beta, evec, ldevec)
          integer    n, neval, lda, ldb, ldevec
          double precision a(lda,*), b(ldb,*), beta(*)
          complex    *16 alpha(*), evec(ldevec,*)
        end function
      end interface

      interface
        double precision function dg2irg (n, neval, a, lda,b, ldb, alpha, beta, evec, ldevec, wk)
          integer    n, neval, lda, ldb, ldevec
          double precision a(lda,*), b(ldb,*), beta(*)
          complex    *16 alpha(*), evec(ldevec,*), wk(n,*)
        end function
      end interface

      interface
        subroutine dgvlcg (n, a, lda, b, ldb, alpha, beta)
          integer    n, lda, ldb
          complex    *16 a(lda,*), b(ldb,*), alpha(*), beta(*)
        end subroutine
      end interface

      interface
        subroutine dg3lcg (n, a, lda, b, ldb, alpha, beta,acopy, bcopy, cwk, wk, iwk)
          integer    n, lda, ldb, iwk(n)
          double precision wk(n)
          complex    *16 a(lda,*), b(ldb,*), alpha(*), beta(*),             &
     &           acopy(n,*), bcopy(n,*), cwk(n)
        end subroutine
      end interface

      interface
        subroutine dgvccg (n, a, lda, b, ldb, alpha, beta,evec, ldevec)
          integer    n, lda, ldb, ldevec
          complex    *16 a(lda,*), b(ldb,*), alpha(*), beta(*),             &
     &           evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine dg6ccg (n, a, lda, b, ldb, alpha, beta,evec, ldevec, acopy, bcopy, cwk, wk, iwk)
          integer    n, lda, ldb, ldevec, iwk(n)
          double precision wk(n)
          complex    *16 a(lda,*), b(ldb,*), alpha(*), beta(*),             &
     &           evec(ldevec,*), acopy(n,*), bcopy(n,*), cwk(n)
        end subroutine
      end interface

      interface
        double precision function dgpicg (n, neval, a, lda,b, ldb, alpha, beta, evec, ldevec)
          integer    n, neval, lda, ldb, ldevec
          complex    *16 a(lda,*), b(ldb,*), alpha(*), beta(*),             &
     &           evec(ldevec,*)
        end function
      end interface

      interface
        double precision function dg2icg (n, neval, a, lda,b, ldb, alpha, beta, evec, ldevec, wk)
          integer    n, neval, lda, ldb, ldevec
          complex    *16 a(lda,*), b(ldb,*), alpha(*), beta(*),             &
     &           evec(ldevec,*), wk(n,*)
        end function
      end interface

      interface
        subroutine dgvcsp (n, a, lda, b, ldb, eval, evec,ldevec)
          integer    n, lda, ldb, ldevec
          double precision a(lda,*), b(ldb,*), eval(*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine dg3csp (n, a, lda, b, ldb, eval, evec,ldevec, iwk, wk1, wk2)
          integer    n, lda, ldb, ldevec, iwk(*)
          double precision a(lda,*), b(ldb,*), eval(*), evec(ldevec,*),     &
     &           wk1(n,*), wk2(n+1,*)
        end subroutine
      end interface

      interface
        double precision function dgpisp (n, neval, a, lda,b, ldb, eval, evec, ldevec)
          integer    n, neval, lda, ldb, ldevec
          double precision a(lda,*), b(ldb,*), eval(*), evec(ldevec,*)
        end function
      end interface

      interface
        double precision function dg2isp (n, neval, a, lda,b, ldb, eval, evec, ldevec, work)
          integer    n, neval, lda, ldb, ldevec
          double precision a(lda,*), b(ldb,*), eval(*), evec(ldevec,*),     &
     &           work(n,*)
        end function
      end interface

!
!     Chapter 3:  Interpolation and Approximation
!

      interface
        subroutine dcsiez (ndata, xdata, fdata, n, xvec,value)
          integer    ndata, n
          double precision xdata(*), fdata(*), xvec(*), value(*)
        end subroutine
      end interface

      interface
        subroutine dc2iez (ndata, xdata, fdata, n, xvec,value, iwk, wk1, wk2)
          integer    ndata, n, iwk(*)
          double precision xdata(*), fdata(*), xvec(*), value(*), wk1(*),   &
     &           wk2(*)
        end subroutine
      end interface

      interface
        subroutine dcsint (ndata, xdata, fdata, break,cscoef)
          integer    ndata
          double precision xdata(*), fdata(*), break(*), cscoef(4,*)
        end subroutine
      end interface

      interface
        subroutine dc2int (ndata, xdata, fdata, break,cscoef, ipvt)
          integer    ndata, ipvt(*)
          double precision xdata(*), fdata(*), break(*), cscoef(4,*)
        end subroutine
      end interface

      interface
        subroutine dcsdec (ndata, xdata, fdata, ileft, dleft,iright, dright, break, cscoef)
          integer    ndata, ileft, iright
          double precision dleft, dright, xdata(*), fdata(*), break(*),     &
     &           cscoef(4,*)
        end subroutine
      end interface

      interface
        subroutine dc2dec (ndata, xdata, fdata, ileft, dleft,iright, dright, break, cscoef, ipvt)
          integer    ndata, ileft, iright, ipvt(*)
          double precision dleft, dright, xdata(*), fdata(*), break(*),     &
     &           cscoef(4,*)
        end subroutine
      end interface

      interface
        subroutine dcsher (ndata, xdata, fdata, dfdata,break, cscoef)
          integer    ndata
          double precision xdata(*), fdata(*), dfdata(*), break(*),         &
     &           cscoef(4,*)
        end subroutine
      end interface

      interface
        subroutine dc2her (ndata, xdata, fdata, dfdata,break, cscoef, ipvt)
          integer    ndata, ipvt(*)
          double precision xdata(*), fdata(*), dfdata(*), break(*),         &
     &           cscoef(4,*)
        end subroutine
      end interface

      interface
        subroutine dcsakm (ndata, xdata, fdata, break,cscoef)
          integer    ndata
          double precision xdata(*), fdata(*), break(*), cscoef(4,*)
        end subroutine
      end interface

      interface
        subroutine dc2akm (ndata, xdata, fdata, break,cscoef, ipvt)
          integer    ndata, ipvt(*)
          double precision xdata(*), fdata(*), break(*), cscoef(4,*)
        end subroutine
      end interface

      interface
        subroutine dcscon (ndata, xdata, fdata, ibreak,break, cscoef)
          integer    ndata, ibreak
          double precision xdata(*), fdata(*), break(*), cscoef(4,*)
        end subroutine
      end interface

      interface
        subroutine dc2con (ndata, xdata, fdata, ibreak,break, cscoef, itmax, xsrt, fsrt, a, y, divd, id, wk)
          integer    ndata, ibreak, itmax, id(*)
          double precision xdata(*), fdata(*), break(*), cscoef(4,*),       &
     &           xsrt(*), fsrt(*), a(*), y(*), divd(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dcsper (ndata, xdata, fdata, break,cscoef)
          integer    ndata
          double precision xdata(*), fdata(*), break(*), cscoef(4,*)
        end subroutine
      end interface

      interface
        subroutine dc2per (ndata, xdata, fdata, break,cscoef, work, ipvt)
          integer    ndata, ipvt(*)
          double precision xdata(*), fdata(*), break(*), cscoef(4,*),       &
     &           work(6,*)
        end subroutine
      end interface

      interface
        double precision function dcsval (x, nintv, break,cscoef)
          integer    nintv
          double precision x, break(*), cscoef(4,*)
        end function
      end interface

      interface
        double precision function dcsder (ideriv, x, nintv,break, cscoef)
          integer    ideriv, nintv
          double precision x, break(*), cscoef(4,*)
        end function
      end interface

      interface
        subroutine dcs1gd (ideriv, n, xvec, nintv, break,cscoef, value)
          integer    ideriv, n, nintv
          double precision xvec(*), break(*), cscoef(4,*), value(*)
        end subroutine
      end interface

      interface
        subroutine dc21gd (ideriv, n, xvec, nintv, break,cscoef, value, left, h, p)
          integer    ideriv, n, nintv, left(*)
          double precision xvec(*), break(*), cscoef(4,*), value(*), h(*),  &
     &           p(*)
        end subroutine
      end interface

      interface
        double precision function dcsitg (a, b, nintv, break,cscoef)
          integer    nintv
          double precision a, b, break(*), cscoef(4,*)
        end function
      end interface

      interface
        subroutine dsplez (ndata, xdata, fdata, itype, ider,n, xvec, value)
          integer    ndata, itype, ider, n
          double precision xdata(*), fdata(*), xvec(*), value(*)
        end subroutine
      end interface

      interface
        subroutine ds2lez (ndata, xdata, fdata, itype, ider,n, xvec, value, wk, iwk)
          integer    ndata, itype, ider, n, iwk(*)
          double precision xdata(*), fdata(*), xvec(*), value(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dbsint (ndata, xdata, fdata, korder,xknot, bscoef)
          integer    ndata, korder
          double precision xdata(*), fdata(*), xknot(*), bscoef(*)
        end subroutine
      end interface

      interface
        subroutine db2int (ndata, xdata, fdata, korder,xknot, bscoef, work1, work2, work3, iwork)
          integer    ndata, korder, iwork(*)
          double precision xdata(*), fdata(*), xknot(*), bscoef(*),         &
     &           work1(*), work2(*), work3(*)
        end subroutine
      end interface

      interface
        subroutine dbsnak (ndata, xdata, korder, xknot)
          integer    ndata, korder
          double precision xdata(*), xknot(*)
        end subroutine
      end interface

      interface
        subroutine db2nak (ndata, xdata, korder, xknot, xsrt,iwk)
          integer    ndata, korder, iwk(*)
          double precision xdata(*), xknot(*), xsrt(*)
        end subroutine
      end interface

      interface
        subroutine dbsopk (ndata, xdata, korder, xknot)
          integer    ndata, korder
          double precision xdata(*), xknot(*)
        end subroutine
      end interface

      interface
        subroutine db2opk (ndata, xdata, korder, xknot,maxit, wk, iwk)
          integer    ndata, korder, maxit, iwk(*)
          double precision xdata(*), xknot(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dbs2in (nxdata, xdata, nydata, ydata,fdata, ldf, kxord, kyord, xknot, yknot, bscoef)
          integer    nxdata, nydata, ldf, kxord, kyord
          double precision xdata(*), ydata(*), fdata(*), xknot(*),          &
     &           yknot(*), bscoef(*)
        end subroutine
      end interface

      interface
        subroutine db22in (nxdata, xdata, nydata, ydata,fdata, ldf, kxord, kyord, xknot, yknot, bscoef, wk,iwk)
          integer    nxdata, nydata, ldf, kxord, kyord, iwk(*)
          double precision xdata(*), ydata(*), fdata(ldf,*), xknot(*),      &
     &           yknot(*), bscoef(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dbs3in (nxdata, xdata, nydata, ydata,nzdata, zdata, fdata, ldf, mdf, kxord, kyord, kzord,xknot, yknot, zknot, bscoef)
          integer    nxdata, nydata, nzdata, ldf, mdf, kxord, kyord, kzord
          double precision xdata(*), ydata(*), zdata(*), fdata(*),          &
     &           xknot(*), yknot(*), zknot(*), bscoef(*)
        end subroutine
      end interface

      interface
        subroutine db23in (nxdata, xdata, nydata, ydata,nzdata, zdata, fdata, ldf, mdf, kxord, kyord, kzord,xknot, yknot, zknot, bscoef, wk, iwk)
          integer    nxdata, nydata, nzdata, ldf, mdf, kxord, kyord,        &
     &           kzord, iwk(*)
          double precision xdata(*), ydata(*), zdata(*), fdata(ldf,mdf,*),  &
     &           xknot(*), yknot(*), zknot(*), bscoef(*), wk(*)
        end subroutine
      end interface

      interface
        double precision function dbsval (x, korder, xknot,ncoef, bscoef)
          integer    korder, ncoef
          double precision x, xknot(*), bscoef(*)
        end function
      end interface

      interface
        double precision function db2val (x, korder, xknot,ncoef, bscoef, aj, dl, dr)
          integer    korder, ncoef
          double precision x, xknot(*), bscoef(*), aj(*), dl(*), dr(*)
        end function
      end interface

      interface
        double precision function dbsder (ideriv, x, korder,xknot, ncoef, bscoef)
          integer    ideriv, korder, ncoef
          double precision x, xknot(*), bscoef(*)
        end function
      end interface

      interface
        double precision function db2der (ideriv, x, korder,xknot, ncoef, bscoef, aj, dl, dr)
          integer    ideriv, korder, ncoef
          double precision x, xknot(*), bscoef(*), aj(*), dl(*), dr(*)
        end function
      end interface

      interface
        subroutine dbs1gd (ideriv, n, xvec, korder, xknot,ncoef, bscoef, value)
          integer    ideriv, n, korder, ncoef
          double precision xvec(*), xknot(*), bscoef(*), value(*)
        end subroutine
      end interface

      interface
        subroutine db21gd (ideriv, n, xvec, korder, xknot,ncoef, bscoef, value, ppcoef, break, left, h, p, wk)
          integer    ideriv, n, korder, ncoef, left(*)
          double precision xvec(*), xknot(*), bscoef(*), value(*),          &
     &           ppcoef(korder,*), break(*), h(*), p(*), wk(*)
        end subroutine
      end interface

      interface
        double precision function dbsitg (a, b, korder,xknot, ncoef, bscoef)
          integer    korder, ncoef
          double precision a, b, xknot(*), bscoef(*)
        end function
      end interface

      interface
        double precision function db2itg (a, b, korder,xknot, ncoef, bscoef, tcoef, aj, dl, dr)
          integer    korder, ncoef
          double precision a, b, xknot(*), bscoef(*), tcoef(0:*), aj(*),    &
     &           dl(*), dr(*)
        end function
      end interface

      interface
        double precision function dbs2vl (x, y, kxord, kyord,xknot, yknot, nxcoef, nycoef, bscoef)
          integer    kxord, kyord, nxcoef, nycoef
          double precision x, y, xknot(*), yknot(*), bscoef(*)
        end function
      end interface

      interface
        double precision function db22vl (x, y, kxord, kyord,xknot, yknot, nxcoef, nycoef, bscoef, wk)
          integer    kxord, kyord, nxcoef, nycoef
          double precision x, y, xknot(*), yknot(*), bscoef(*), wk(*)
        end function
      end interface

      interface
        double precision function dbs2dr (ixder, iyder, x, y,kxord, kyord, xknot, yknot, nxcoef, nycoef, bscoef)
          integer    ixder, iyder, kxord, kyord, nxcoef, nycoef
          double precision x, y, xknot(*), yknot(*), bscoef(*)
        end function
      end interface

      interface
        double precision function db22dr (ixder, iyder, x, y,kxord, kyord, xknot, yknot, nxcoef, nycoef, bscoef, wk)
          integer    ixder, iyder, kxord, kyord, nxcoef, nycoef
          double precision x, y, xknot(*), yknot(*), bscoef(nxcoef,*),      &
     &           wk(*)
        end function
      end interface

      interface
        subroutine dbs2gd (ixder, iyder, nx, xvec, ny, yvec,kxord, kyord, xknot, yknot, nxcoef, nycoef, bscoef,value, ldvalu)
          integer    ixder, iyder, nx, ny, kxord, kyord, nxcoef, nycoef,    &
     &           ldvalu
          double precision xvec(*), yvec(*), xknot(*), yknot(*),            &
     &           bscoef(*), value(ldvalu,*)
        end subroutine
      end interface

      interface
        subroutine db22gd (ixder, iyder, nx, xvec, ny, yvec,kxord, kyord, xknot, yknot, nxcoef, nycoef, bscoef,value, ldvalu, leftx, lefty, a, b, dbiatx, dbiaty, bx,by)
          integer    ixder, iyder, nx, ny, kxord, kyord, nxcoef, nycoef,    &
     &           ldvalu, leftx(*), lefty(*)
          double precision xvec(*), yvec(*), xknot(*), yknot(*),            &
     &           bscoef(nxcoef,*), value(ldvalu,*), a(kxord,*),         &
     &           b(kyord,*), dbiatx(kxord,*), dbiaty(kyord,*),          &
     &           bx(kxord,*), by(kyord,*)
        end subroutine
      end interface

      interface
        double precision function dbs2ig (a, b, c, d, kxord,kyord, xknot, yknot, nxcoef, nycoef, bscoef)
          integer    kxord, kyord, nxcoef, nycoef
          double precision a, b, c, d, xknot(*), yknot(*), bscoef(*)
        end function
      end interface

      interface
        double precision function db22ig (a, b, c, d, kxord,kyord, xknot, yknot, nxcoef, nycoef, bscoef, wk)
          integer    kxord, kyord, nxcoef, nycoef
          double precision a, b, c, d, xknot(*), yknot(*),                  &
     &           bscoef(nxcoef,*), wk(*)
        end function
      end interface

      interface
        double precision function dbs3vl (x, y, z, kxord,kyord, kzord, xknot, yknot, zknot, nxcoef, nycoef,nzcoef, bscoef)
          integer    kxord, kyord, kzord, nxcoef, nycoef, nzcoef
          double precision x, y, z, xknot(*), yknot(*), zknot(*), bscoef(*)
        end function
      end interface

      interface
        double precision function db23vl (x, y, z, kxord,kyord, kzord, xknot, yknot, zknot, nxcoef, nycoef,nzcoef, bscoef, wk)
          integer    kxord, kyord, kzord, nxcoef, nycoef, nzcoef
          double precision x, y, z, xknot(*), yknot(*), zknot(*),           &
     &           bscoef(nxcoef,nycoef,*), wk(*)
        end function
      end interface

      interface
        double precision function dbs3dr (ixder, iyder,izder, x, y, z, kxord, kyord, kzord, xknot, yknot,zknot, nxcoef, nycoef, nzcoef, bscoef)
          integer    ixder, iyder, izder, kxord, kyord, kzord, nxcoef,      &
     &           nycoef, nzcoef
          double precision x, y, z, xknot(*), yknot(*), zknot(*), bscoef(*)
        end function
      end interface

      interface
        double precision function db23dr (ixder, iyder,izder, x, y, z, kxord, kyord, kzord, xknot, yknot,zknot, nxcoef, nycoef, nzcoef, bscoef, wk)
          integer    ixder, iyder, izder, kxord, kyord, kzord, nxcoef,      &
     &           nycoef, nzcoef
          double precision x, y, z, xknot(*), yknot(*), zknot(*),           &
     &           bscoef(nxcoef,nycoef,*), wk(*)
        end function
      end interface

      interface
        subroutine dbs3gd (ixder, iyder, izder, nx, xvec, ny,yvec, nz, zvec, kxord, kyord, kzord, xknot, yknot,zknot, nxcoef, nycoef, nzcoef, bscoef, value, ldvalu,mdvalu)
          integer    ixder, iyder, izder, nx, ny, nz, kxord, kyord, kzord,  &
     &           nxcoef, nycoef, nzcoef, ldvalu, mdvalu
          double precision xvec(*), yvec(*), zvec(*), xknot(*), yknot(*),   &
     &           zknot(*), bscoef(*), value(ldvalu,mdvalu,*)
        end subroutine
      end interface

      interface
        subroutine db23gd (ixder, iyder, izder, nx, xvec, ny,yvec, nz, zvec, kxord, kyord, kzord, xknot, yknot,zknot, nxcoef, nycoef, nzcoef, bscoef, value, ldvalu,mdvalu, leftx, lefty, leftz, a, b, c, dbiatx, dbiaty,dbiatz, bx, by, bz)
          integer    ixder, iyder, izder, nx, ny, nz, kxord, kyord, kzord,  &
     &           nxcoef, nycoef, nzcoef, ldvalu, mdvalu, leftx(*),      &
     &           lefty(*), leftz(*)
          double precision xvec(*), yvec(*), zvec(*), xknot(*), yknot(*),   &
     &           zknot(*), bscoef(nxcoef,nycoef,*),                     &
     &           value(ldvalu,mdvalu,*), a(kxord,*), b(kyord,*),        &
     &           c(kzord,*), dbiatx(kxord,*), dbiaty(kyord,*),          &
     &           dbiatz(kzord,*), bx(kxord,*), by(kyord,*), bz(kzord,*)
        end subroutine
      end interface

      interface
        double precision function dbs3ig (a, b, c, d, e, f,kxord, kyord, kzord, xknot, yknot, zknot, nxcoef,nycoef, nzcoef, bscoef)
          integer    kxord, kyord, kzord, nxcoef, nycoef, nzcoef
          double precision a, b, c, d, e, f, xknot(*), yknot(*), zknot(*),  &
     &           bscoef(*)
        end function
      end interface

      interface
        double precision function db23ig (a, b, c, d, e, f,kxord, kyord, kzord, xknot, yknot, zknot, nxcoef,nycoef, nzcoef, bscoef, wk)
          integer    kxord, kyord, kzord, nxcoef, nycoef, nzcoef
          double precision a, b, c, d, e, f, xknot(*), yknot(*), zknot(*),  &
     &           bscoef(*), wk(*)
        end function
      end interface

      interface
        subroutine dbscpp (korder, xknot, ncoef, bscoef,nppcf, break, ppcoef)
          integer    korder, ncoef, nppcf
          double precision xknot(*), bscoef(*), break(*), ppcoef(*)
        end subroutine
      end interface

      interface
        subroutine db2cpp (korder, xknot, ncoef, bscoef,nppcf, break, ppcoef, wk)
          integer    korder, ncoef, nppcf
          double precision xknot(*), bscoef(*), break(*), ppcoef(*), wk(*)
        end subroutine
      end interface

      interface
        double precision function dppval (x, korder, nintv,break, ppcoef)
          integer    korder, nintv
          double precision x, break(*), ppcoef(*)
        end function
      end interface

      interface
        double precision function dppder (ideriv, x, korder,nintv, break, ppcoef)
          integer    ideriv, korder, nintv
          double precision x, break(*), ppcoef(korder,*)
        end function
      end interface

      interface
        subroutine dpp1gd (ideriv, n, xvec, korder, nintv,break, ppcoef, value)
          integer    ideriv, n, korder, nintv
          double precision xvec(*), break(*), ppcoef(*), value(*)
        end subroutine
      end interface

      interface
        subroutine dp21gd (ideriv, n, xvec, korder, nintv,break, ppcoef, value, left, h, p)
          integer    ideriv, n, korder, nintv, left(*)
          double precision xvec(*), break(*), ppcoef(korder,*), value(*),   &
     &           h(*), p(*)
        end subroutine
      end interface

      interface
        double precision function dppitg (a, b, korder,nintv, break, ppcoef)
          integer    korder, nintv
          double precision a, b, break(*), ppcoef(korder,*)
        end function
      end interface

      interface
        double precision function dqdval (x, ndata, xdata,fdata, check)
          integer    ndata
          double precision x, xdata(*), fdata(*)
          logical    check
        end function
      end interface

      interface
        double precision function dqdder (ideriv, x, ndata,xdata, fdata, check)
          integer    ideriv, ndata
          double precision x, xdata(*), fdata(*)
          logical    check
        end function
      end interface

      interface
        double precision function dqd2vl (x, y, nxdata,xdata, nydata, ydata, fdata, ldf, check)
          integer    nxdata, nydata, ldf
          double precision x, y, xdata(*), ydata(*), fdata(*)
          logical    check
        end function
      end interface

      interface
        double precision function dqd2dr (ixder, iyder, x, y,nxdata, xdata, nydata, ydata, fdata, ldf, check)
          integer    ixder, iyder, nxdata, nydata, ldf
          double precision x, y, xdata(*), ydata(*), fdata(ldf,*)
          logical    check
        end function
      end interface

      interface
        double precision function dqd3vl (x, y, z, nxdata,xdata, nydata, ydata, nzdata, zdata, fdata, ldf, mdf,check)
          integer    nxdata, nydata, nzdata, ldf, mdf
          double precision x, y, z, xdata(*), ydata(*), zdata(*), fdata(*)
          logical    check
        end function
      end interface

      interface
        double precision function dqd3dr (ixder, iyder,izder, x, y, z, nxdata, xdata, nydata, ydata, nzdata,zdata, fdata, ldf, mdf, check)
          integer    ixder, iyder, izder, nxdata, nydata, nzdata, ldf, mdf
          double precision x, y, z, xdata(*), ydata(*), zdata(*),           &
     &           fdata(ldf,mdf,*)
          logical    check
        end function
      end interface

      interface
        subroutine dsurf (ndata, xydata, fdata, nxout, nyout,xout, yout, sur, ldsur)
          integer    ndata, nxout, nyout, ldsur
          double precision xydata(2,*), fdata(*), xout(*), yout(*), sur(*)
        end subroutine
      end interface

      interface
        subroutine ds2rf (ndata, xydata, fdata, nxout, nyout,xout, yout, sur, ldsur, iwk, wk)
          integer    ndata, nxout, nyout, ldsur, iwk(*)
          double precision xydata(2,*), fdata(*), xout(*), yout(*),         &
     &           sur(ldsur,*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dfnlsq (f, intcep, nbasis, ndata, xdata,fdata, iwt, weight, a, sse)
          integer    intcep, nbasis, ndata, iwt
          double precision f, sse, xdata(*), fdata(*), weight(*), a(*)
          external   f
        end subroutine
      end interface

      interface
        subroutine df2lsq (f, intcep, nbasis, ndata, xdata,fdata, iwt, weight, a, sse, wk)
          integer    intcep, nbasis, ndata, iwt
          double precision f, sse, xdata(*), fdata(*), weight(*), a(*),     &
     &           wk(*)
          external   f
        end subroutine
      end interface

      interface
        subroutine dbslsq (ndata, xdata, fdata, weight,korder, xknot, ncoef, bscoef)
          integer    ndata, korder, ncoef
          double precision xdata(*), fdata(*), weight(*), xknot(*),         &
     &           bscoef(*)
        end subroutine
      end interface

      interface
        subroutine db2lsq (ndata, xdata, fdata, weight,korder, xknot, ncoef, bscoef, wk, xsort, fsort, wsort,ipvt)
          integer    ndata, korder, ncoef, ipvt(*)
          double precision xdata(*), fdata(*), weight(*), xknot(*),         &
     &           bscoef(*), wk(*), xsort(*), fsort(*), wsort(*)
        end subroutine
      end interface

      interface
        subroutine dbsvls (ndata, xdata, fdata, weight,korder, ncoef, xguess, xknot, bscoef, ssq)
          integer    ndata, korder, ncoef
          double precision ssq, xdata(*), fdata(*), weight(*), xguess(*),   &
     &           xknot(*), bscoef(*)
        end subroutine
      end interface

      interface
        subroutine db2vls (ndata, xdata, fdata, weight,korder, ncoef, xguess, xknot, bscoef, ssq, iwk, wk)
          integer    ndata, korder, ncoef, iwk(*)
          double precision ssq, xdata(*), fdata(*), weight(*), xguess(*),   &
     &           xknot(*), bscoef(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dconft (ndata, xdata, fdata, weight,nxval, xval, nhard, ider, itype, bl, bu, korder, xknot,ncoef, bscoef)
          integer    ndata, nxval, nhard, korder, ncoef, ider(*), itype(*)
          double precision xdata(*), fdata(*), weight(*), xval(*), bl(*),   &
     &           bu(*), xknot(*), bscoef(*)
        end subroutine
      end interface

      interface
        subroutine dc2nft (ndata, xdata, fdata, weight,nxval, xval, nhard, ider, itype, bl, bu, korder, xknot,ncoef, bscoef, h, g, a, rhs, wk, iperm, iwk)
          integer    ndata, nhard, nxval, korder, ncoef, ider(*),           &
     &           itype(*), iperm(*), iwk(*)
          double precision xdata(*), fdata(*), weight(*), xval(*), bl(*),   &
     &           bu(*), xknot(*), bscoef(*), h(*), g(*), a(*), rhs(*),  &
     &           wk(*)
        end subroutine
      end interface

      interface
        subroutine dbsls2 (nxdata, xdata, nydata, ydata,fdata, ldf, kxord, kyord, xknot, yknot, nxcoef, nycoef,xweigh, yweigh, bscoef)
          integer    nxdata, nydata, ldf, kxord, kyord, nxcoef, nycoef
          double precision xdata(*), ydata(*), fdata(*), xknot(*),          &
     &           yknot(*), xweigh(*), yweigh(*), bscoef(*)
        end subroutine
      end interface

      interface
        subroutine db2ls2 (nxdata, xdata, nydata, ydata,fdata, ldf, kxord, kyord, xknot, yknot, nxcoef, nycoef,xweigh, yweigh, bscoef, wk)
          integer    nxdata, nydata, ldf, kxord, kyord, nxcoef, nycoef
          double precision xdata(*), ydata(*), fdata(*), xknot(*),          &
     &           yknot(*), xweigh(*), yweigh(*), bscoef(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dbsls3 (nxdata, xdata, nydata, ydata,nzdata, zdata, fdata, ldfdat, mdfdat, kxord, kyord,kzord, xknot, yknot, zknot, nxcoef, nycoef, nzcoef,xweigh, yweigh, zweigh, bscoef)
          integer    nxdata, nydata, nzdata, ldfdat, mdfdat, kxord, kyord,  &
     &           kzord, nxcoef, nycoef, nzcoef
          double precision xdata(*), ydata(*), zdata(*), fdata(*),          &
     &           xknot(*), yknot(*), zknot(*), xweigh(*), yweigh(*),    &
     &           zweigh(*), bscoef(nxcoef,nycoef,*)
        end subroutine
      end interface

      interface
        subroutine db2ls3 (nxdata, xdata, nydata, ydata,nzdata, zdata, fdata, ldfdat, mdfdat, kxord, kyord,kzord, xknot, yknot, zknot, nxcoef, nycoef, nzcoef,xweigh, yweigh, zweigh, bscoef, wk)
          integer    nxdata, nydata, nzdata, ldfdat, mdfdat, kxord, kyord,  &
     &           kzord, nxcoef, nycoef, nzcoef
          double precision zknot(*), xdata(*), ydata(*), zdata(*),          &
     &           fdata(*), xknot(*), yknot(*), xweigh(*), yweigh(*),    &
     &           zweigh(*), bscoef(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dcssed (ndata, xdata, fdata, dis, sc,maxit, sdata)
          integer    ndata, maxit
          double precision dis, sc, xdata(*), fdata(*), sdata(*)
        end subroutine
      end interface

      interface
        subroutine dc2sed (ndata, xdata, fdata, dis, sc,maxit, sdata, wk, iwk)
          integer    ndata, maxit, iwk(*)
          double precision dis, sc, xdata(*), fdata(*), sdata(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dcssmh (ndata, xdata, fdata, weight,smpar, break, cscoef)
          integer    ndata
          double precision smpar, xdata(*), fdata(*), weight(*), break(*),  &
     &           cscoef(4,*)
        end subroutine
      end interface

      interface
        subroutine dc2smh (ndata, xdata, fdata, weight,smpar, break, cscoef, wk, iwk)
          integer    ndata, iwk(*)
          double precision smpar, xdata(*), fdata(*), weight(*), break(*),  &
     &           cscoef(4,*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dcsscv (ndata, xdata, fdata, iequal,break, cscoef)
          integer    ndata, iequal
          double precision xdata(*), fdata(*), break(*), cscoef(4,*)
        end subroutine
      end interface

      interface
        subroutine dc2scv (ndata, xdata, fdata, iequal,break, cscoef, wk, ywk, iwk)
          integer    ndata, iequal, iwk(*)
          double precision xdata(*), fdata(*), break(*), cscoef(4,*),       &
     &           wk(0:ndata+1,7), ywk(*)
        end subroutine
      end interface

      interface
        subroutine dratch (f, phi, weight, a, b, n, m, p, q,error)
          integer    n, m
          double precision f, phi, weight, a, b, error, p(*), q(*)
          external   f, phi, weight
        end subroutine
      end interface

      interface
        subroutine dr2tch (f, phi, weight, a, b, n, m, p, q,error, itmax, iwk, wk)
          integer    n, m, itmax, iwk(*)
          double precision f, phi, weight, a, b, error, p(*), q(*), wk(*)
          external   f, phi, weight
        end subroutine
      end interface

!
!     Chapter 4:  Integration and Differentiation
!

      interface
        subroutine dqdag (f, a, b, errabs, errrel, irule,result, errest)
          integer    irule
          double precision f, a, b, errabs, errrel, result, errest
          external   f
        end subroutine
      end interface

      interface
        subroutine dq2ag (f, a, b, errabs, errrel, irule,result, errest, maxsub, neval, nsubin, alist, blist,rlist, elist, iord)
          integer    irule, maxsub, neval, nsubin, iord(*)
          double precision f, a, b, errabs, errrel, result, errest,         &
     &           alist(*), blist(*), rlist(*), elist(*)
          external   f
        end subroutine
      end interface

      interface
        subroutine dqdagp (f, a, b, npts, points, errabs,errrel, result, errest)
          integer    npts
          double precision f, a, b, errabs, errrel, result, errest,         &
     &           points(*)
          external   f
        end subroutine
      end interface

      interface
        subroutine dq2agp (f, a, b, npts, points, errabs,errrel, result, errest, maxsub, neval, nsubin, alist,blist, rlist, elist, iord, level, wk, iwk)
          integer    npts, maxsub, neval, nsubin, iord(*), level(*), iwk(*)
          double precision f, a, b, errabs, errrel, result, errest,         &
     &           points(*), alist(*), blist(*), rlist(*), elist(*),     &
     &           wk(*)
          external   f
        end subroutine
      end interface

      interface
        subroutine dqdagi (f, bound, interv, errabs, errrel,result, errest)
          integer    interv
          double precision f, bound, errabs, errrel, result, errest
          external   f
        end subroutine
      end interface

      interface
        subroutine dq2agi (f, bound, interv, errabs, errrel,result, errest, maxsub, neval, nsubin, alist, blist,rlist, elist, iord)
          integer    interv, maxsub, neval, nsubin, iord(*)
          double precision f, bound, errabs, errrel, result, errest,        &
     &           alist(*), blist(*), rlist(*), elist(*)
        end subroutine
      end interface

      interface
        subroutine dqdawo (f, a, b, iweigh, omega, errabs,errrel, result, errest)
          integer    iweigh
          double precision f, a, b, omega, errabs, errrel, result, errest
          external   f
        end subroutine
      end interface

      interface
        subroutine dq2awo (f, a, b, iweigh, omega, errabs,errrel, result, errest, maxsub, maxcby, neval, nsubin,alist, blist, rlist, elist, iord, nnlog, wk)
          integer    iweigh, maxsub, maxcby, neval, nsubin, iord(*),        &
     &           nnlog(*)
          double precision f, a, b, omega, errabs, errrel, result, errest,  &
     &           alist(*), blist(*), rlist(*), elist(*), wk(*)
          external   f
        end subroutine
      end interface

      interface
        subroutine dqdawf (f, a, iweigh, omega, errabs,result, errest)
          integer    iweigh
          double precision f, a, omega, errabs, result, errest
          external   f
        end subroutine
      end interface

      interface
        subroutine dq2awf (f, a, iweigh, omega, errabs,result, errest, maxcyl, maxsub, maxcby, neval, ncycle,rslist, erlist, ierlst, nsubin, wk, iwk)
          integer    iweigh, maxcyl, maxsub, maxcby, neval, ncycle,         &
     &           nsubin, ierlst(*), iwk(2,*)
          double precision f, a, omega, errabs, result, errest, rslist(*),  &
     &           erlist(*), wk(*)
          external   f
        end subroutine
      end interface

      interface
        subroutine dqdaws (f, a, b, iweigh, alpha, beta,errabs, errrel, result, errest)
          integer    iweigh
          double precision f, a, b, alpha, beta, errabs, errrel, result,    &
     &           errest
          external   f
        end subroutine
      end interface

      interface
        subroutine dq2aws (f, a, b, iweigh, alpha, beta,errabs, errrel, result, errest, maxsub, neval, nsubin,alist, blist, rlist, elist, iord)
          integer    iweigh, maxsub, neval, nsubin, iord(*)
          double precision f, a, b, alpha, beta, errabs, errrel, result,    &
     &           errest, alist(*), blist(*), rlist(*), elist(*)
          external   f
        end subroutine
      end interface

      interface
        subroutine dqdawc (f, a, b, c, errabs, errrel,result, errest)
          double precision f, a, b, c, errabs, errrel, result, errest
          external   f
        end subroutine
      end interface

      interface
        subroutine dq2awc (f, a, b, c, errabs, errrel,result, errest, maxsub, neval, nsubin, alist, blist,rlist, elist, iord)
          integer    maxsub, neval, nsubin, iord(*)
          double precision f, a, b, c, errabs, errrel, result, errest,      &
     &           alist(*), blist(*), rlist(*), elist(*)
          external   f
        end subroutine
      end interface

      interface
        subroutine dqdng (f, a, b, errabs, errrel, result,errest)
          double precision f, a, b, errabs, errrel, result, errest
          external   f
        end subroutine
      end interface

      interface
        subroutine dtwodq (f, a, b, g, h, errabs, errrel,irule, result, errest)
          integer    irule
          double precision f, a, b, g, h, errabs, errrel, result, errest
          external   f, g, h
        end subroutine
      end interface

      interface
        subroutine dt2odq (user, a, b, g, h, errabs, errrel,irule, result, errest, maxsub, neval, nsubin, alist,blist, rlist, elist, iord, wk, iwk)
          integer    irule, maxsub, neval, nsubin, iord(*), iwk(*)
          double precision user, a, b, g, h, errabs, errrel, result,        &
     &           errest, alist(*), blist(*), rlist(*), elist(*), wk(*)
          external   user, g, h
        end subroutine
      end interface

      interface
        subroutine dqand (f, n, a, b, errabs, errrel, maxfcn,result, errest)
          integer    n, maxfcn
          double precision f, errabs, errrel, result, errest, a(*), b(*)
          external   f
        end subroutine
      end interface

      interface
        subroutine dgqrul (n, iweigh, alpha, beta, nfix,qxfix, qx, qw)
          integer    n, iweigh, nfix
          double precision alpha, beta, qxfix(*), qx(*), qw(*)
        end subroutine
      end interface

      interface
        subroutine dg2rul (n, iweigh, alpha, beta, nfix,qxfix, qx, qw, wk)
          integer    n, iweigh, nfix
          double precision alpha, beta, qxfix(*), qx(*), qw(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dgqrcf (n, b, c, nfix, qxfix, qx, qw)
          integer    n, nfix
          double precision b(*), c(*), qxfix(*), qx(*), qw(*)
        end subroutine
      end interface

      interface
        subroutine dg2rcf (n, b, c, nfix, qxfix, qx, qw, wk)
          integer    n, nfix
          double precision b(*), c(*), qxfix(*), qx(*), qw(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dreccf (n, iweigh, alpha, beta, b, c)
          integer    n, iweigh
          double precision alpha, beta, b(*), c(*)
        end subroutine
      end interface

      interface
        subroutine drecqr (n, qx, qw, nterm, b, c)
          integer    n, nterm
          double precision qx(*), qw(*), b(*), c(*)
        end subroutine
      end interface

      interface
        subroutine dr2cqr (n, qx, qw, nterm, b, c, wk)
          integer    n, nterm
          double precision qx(*), qw(*), b(*), c(*), wk(2,*)
        end subroutine
      end interface

      interface
        subroutine dfqrul (n, a, b, iweigh, alpha, beta, qx,qw)
          integer    n, iweigh
          double precision a, b, alpha, beta, qx(*), qw(*)
        end subroutine
      end interface

      interface
        subroutine df2rul (n, a, b, iweigh, alpha, beta, qx,qw, wk)
          integer    n, iweigh
          double precision a, b, alpha, beta, qx(*), qw(*), wk(*)
        end subroutine
      end interface

      interface
        double precision function dderiv (fcn, korder, x,bgstep, tol)
          integer    korder
          double precision fcn, x, bgstep, tol
          external   fcn
        end function
      end interface

!
!     Chapter 5:  Differential Equations
!

      interface
        subroutine divprk (ido, neq, fcn, x, xend, tol,param, y)
          integer    ido, neq
          double precision x, xend, tol, param(*), y(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine di2prk (ido, neq, fcn, x, xend, tol,param, y, vnorm, wk)
          integer    ido, neq
          double precision x, xend, tol, param(*), y(*), wk(neq,*)
          external   fcn, vnorm
        end subroutine
      end interface

      interface
        subroutine di3prk (neq, v, y, ymax, enorm)
          integer    neq
          double precision enorm, v(*), y(*), ymax(*)
        end subroutine
      end interface

      interface
        subroutine divpag (ido, neq, fcn, fcnj, a, x, xend,tol, param, y)
          integer    ido, neq
          double precision x, xend, tol, a(*), param(*), y(*)
          external   fcn, fcnj
        end subroutine
      end interface

      interface
        subroutine di2pag (ido, neq, fcn, fcnj, a, x, xend,tol, param, y, ytemp, ymax, error, save1, save2, pw,ipvt, vnorm)
          integer    ido, neq, ipvt(*)
          double precision x, xend, tol, a(*), param(*), y(*), ytemp(*),    &
     &           ymax(*), error(*), save1(*), save2(*), pw(*)
          external   fcn, fcnj, vnorm
        end subroutine
      end interface

      interface
        subroutine dbvpfd (fcneqn, fcnjac, fcnbc, fcnpeq,fcnpbc, neqns, nleft, ncupbc, xleft, xright, pistep,tol, ninit, xinit, yinit, ldyini, linear, print,mxgrid, nfinal, xfinal, yfinal, ldyfin, errest)
          integer    neqns, nleft, ncupbc, ninit, ldyini, mxgrid, nfinal,   &
     &           ldyfin
          double precision xleft, xright, pistep, tol, xinit(*),            &
     &           yinit(ldyini,*), xfinal(*), yfinal(ldyfin,*),          &
     &           errest(*)
          logical    linear, print
          external   fcneqn, fcnjac, fcnbc, fcnpeq, fcnpbc
        end subroutine
      end interface

      interface
        subroutine db2pfd (fcneqn, fcnjac, fcnbc, fcnpeq,fcnpbc, neqns, nleft, ncupbc, xleft, xright, pistep,tol, ninit, xinit, yinit, ldyini, linear, print,mxgrid, nfinal, xfinal, yfinal, ldyfin, errest, rwork,iwork)
          integer    neqns, nleft, ncupbc, ninit, ldyini, mxgrid, nfinal,   &
     &           ldyfin, iwork(*)
          double precision xleft, xright, pistep, tol, xinit(*),            &
     &           yinit(ldyini,*), xfinal(*), yfinal(ldyfin,*),          &
     &           errest(*), rwork(*)
          logical    linear, print
          external   fcneqn, fcnjac, fcnbc, fcnpeq, fcnpbc
        end subroutine
      end interface

      interface
        subroutine dbvpms (fcneqn, fcnjac, fcnbc, neqns,xleft, xright, dtol, btol, maxit, ninit, xinit, yinit,ldyini, nmax, nfinal, xfinal, yfinal, ldyfin)
          integer    neqns, maxit, ninit, ldyini, nmax, nfinal, ldyfin
          double precision xleft, xright, dtol, btol, xinit(*),             &
     &           yinit(ldyini,*), xfinal(*), yfinal(ldyfin,*)
          external   fcneqn, fcnjac, fcnbc
        end subroutine
      end interface

      interface
        subroutine db2pms (fcneqn, fcnjac, fcnbc, neqns,xleft, xright, dtol, btol, maxit, ninit, xinit, yinit,ldyini, nmax, nfinal, xfinal, yfinal, ldyfin, work,iwk)
          integer    neqns, maxit, ninit, ldyini, nmax, nfinal, ldyfin,     &
     &           iwk(*)
          double precision xleft, xright, dtol, btol, xinit(*),             &
     &           yinit(ldyini,*), xfinal(*), yfinal(ldyfin,*), work(*)
          external   fcneqn, fcnjac, fcnbc
        end subroutine
      end interface

      interface
        subroutine ddaspg (neq, t, tout, ido, y, yprime,ddgspg)
          integer    neq, ido
          double precision t, tout, y(*), yprime(*)
          external   ddgspg
        end subroutine
      end interface

      interface
        subroutine dd2spg (neq, t, tout, ido, y, yprime,ddgspg, ddjspg, iwork, rwork)
          integer    neq, ido, iwork(*)
          double precision t, tout, y(*), yprime(*), rwork(*)
          external   ddgspg, ddjspg
        end subroutine
      end interface

      interface
        subroutine ddgspg (n, t, y, ypr, gval)
          integer    n
          double precision t, y(n), ypr(n), gval(n)
        end subroutine
      end interface

      interface
        subroutine ddjspg (n, t, y, ypr, cj, pdg, ldpdg)
          integer    n, ldpdg
          double precision t, cj, y(n), ypr(n), pdg(ldpdg*n)
        end subroutine
      end interface

      interface
        subroutine dmolch (ido, fcnut, fcnbc, npdes, t, tend,nx, xbreak, tol, hinit, y, ldy)
          integer    ido, npdes, nx, ldy
          double precision t, tend, tol, hinit, xbreak(*), y(ldy,*)
          external   fcnut, fcnbc
        end subroutine
      end interface

      interface
        subroutine dm2lch (ido, fcnut, fcnbc, npdes, t, tend,nx, xbreak, tol, hinit, y, ldy, wk, iwk)
          integer    ido, npdes, nx, ldy, iwk(*)
          double precision t, tend, tol, hinit, xbreak(*), y(ldy,*), wk(*)
          external   fcnut, fcnbc
        end subroutine
      end interface

      interface
        subroutine dfps2h (prhs, brhs, coefu, nx, ny, ax, bx,ay, by, ibcty, iorder, u, ldu)
          integer    nx, ny, iorder, ldu, ibcty(*)
          double precision prhs, brhs, coefu, ax, bx, ay, by, u(ldu,*)
          external   prhs, brhs
        end subroutine
      end interface

      interface
        subroutine df2s2h (prhs, brhs, coefu, nx, ny, ax, bx,ay, by, ibcty, iorder, u, ldu, uwork, work)
          integer    nx, ny, iorder, ldu, ibcty(*)
          double precision coefu, ax, bx, ay, by, u(ldu,*), uwork(*),       &
     &           work(*), prhs, brhs
          external   prhs, brhs
        end subroutine
      end interface

      interface
        subroutine dfps3h (prhs, brhs, coefu, nx, ny, nz, ax,bx, ay, by, az, bz, ibcty, iorder, u, ldu, mdu)
          integer    nx, ny, nz, iorder, ldu, mdu, ibcty(*)
          double precision prhs, brhs, coefu, ax, bx, ay, by, az, bz,       &
     &           u(ldu,mdu,*)
          external   prhs, brhs
        end subroutine
      end interface

      interface
        subroutine df2s3h (prhs, brhs, coefu, nx, ny, nz, ax,bx, ay, by, az, bz, ibcty, iorder, u, ldu, mdu, uwork,work)
          integer    nx, ny, nz, iorder, ldu, mdu, ibcty(*)
          double precision prhs, brhs, coefu, ax, bx, ay, by, az, bz,       &
     &           u(ldu,mdu,*), uwork(*), work(*)
          external   prhs, brhs
        end subroutine
      end interface

!
!     Chapter 6:  Transforms
!

      interface
        subroutine dfftcf (n, seq, coef)
          integer    n
          complex    *16 seq(*), coef(*)
        end subroutine
      end interface

      interface
        subroutine df2tcf (n, seq, coef, wfftc, cpy)
          integer    n
          double precision wfftc(*), cpy(*)
          complex    *16 seq(*), coef(*)
        end subroutine
      end interface

      interface
        subroutine dfftcb (n, coef, seq)
          integer    n
          complex    *16 coef(*), seq(*)
        end subroutine
      end interface

      interface
        subroutine df2tcb (n, coef, seq, wfftc, cpy)
          integer    n
          double precision wfftc(*), cpy(*)
          complex    *16 coef(*), seq(*)
        end subroutine
      end interface

      interface
        subroutine dfftci (n, wfftc)
          integer    n
          double precision wfftc(*)
        end subroutine
      end interface

      interface
        subroutine dfsint (n, seq, coef)
          integer    n
          double precision seq(*), coef(*)
        end subroutine
      end interface

      interface
        subroutine df2int (n, seq, coef, wfsin)
          integer    n
          double precision seq(*), coef(*), wfsin(*)
        end subroutine
      end interface

      interface
        subroutine dfsini (n, wfsin)
          integer    n
          double precision wfsin(*)
        end subroutine
      end interface

      interface
        subroutine dfcost (n, seq, coef)
          integer    n
          double precision seq(*), coef(*)
        end subroutine
      end interface

      interface
        subroutine df2ost (n, seq, coef, wfcos)
          integer    n
          double precision seq(*), coef(*), wfcos(*)
        end subroutine
      end interface

      interface
        subroutine dfcosi (n, wfcos)
          integer    n
          double precision wfcos(*)
        end subroutine
      end interface

      interface
        subroutine dqsinf (n, seq, coef)
          integer    n
          double precision seq(*), coef(*)
        end subroutine
      end interface

      interface
        subroutine dq2inf (n, seq, coef, wqsin)
          integer    n
          double precision seq(*), coef(*), wqsin(*)
        end subroutine
      end interface

      interface
        subroutine dqsinb (n, coef, seq)
          integer    n
          double precision coef(*), seq(*)
        end subroutine
      end interface

      interface
        subroutine dq2inb (n, coef, seq, wqsin)
          integer    n
          double precision coef(*), seq(*), wqsin(*)
        end subroutine
      end interface

      interface
        subroutine dqsini (n, wqsin)
          integer    n
          double precision wqsin(*)
        end subroutine
      end interface

      interface
        subroutine dqcosf (n, seq, coef)
          integer    n
          double precision seq(*), coef(*)
        end subroutine
      end interface

      interface
        subroutine dq2osf (n, seq, coef, wqcos)
          integer    n
          double precision seq(*), coef(*), wqcos(*)
        end subroutine
      end interface

      interface
        subroutine dqcosb (n, coef, seq)
          integer    n
          double precision coef(*), seq(*)
        end subroutine
      end interface

      interface
        subroutine dq2osb (n, coef, seq, wqcos)
          integer    n
          double precision coef(*), seq(*), wqcos(*)
        end subroutine
      end interface

      interface
        subroutine dqcosi (n, wqcos)
          integer    n
          double precision wqcos(*)
        end subroutine
      end interface

      interface
        subroutine dfft2d (nra, nca, a, lda, coef, ldcoef)
          integer    nra, nca, lda, ldcoef
          complex    *16 a(lda,*), coef(ldcoef,*)
        end subroutine
      end interface

      interface
        subroutine df2t2d (nra, nca, a, lda, coef, ldcoef,wff1, wff2, cwk, cpy)
          integer    nra, nca, lda, ldcoef
          double precision wff1(*), wff2(*), cpy(*)
          complex    *16 a(lda,*), coef(ldcoef,*), cwk(*)
        end subroutine
      end interface

      interface
        subroutine dfft2b (nrcoef, nccoef, coef, ldcoef, a,lda)
          integer    nrcoef, nccoef, ldcoef, lda
          complex    *16 coef(ldcoef,*), a(lda,*)
        end subroutine
      end interface

      interface
        subroutine df2t2b (nrcoef, nccoef, a, lda, coef,ldcoef, wff1, wff2, cwk, cpy)
          integer    nrcoef, nccoef, lda, ldcoef
          double precision wff1(*), wff2(*), cpy(*)
          complex    *16 a(lda,*), coef(ldcoef,*), cwk(*)
        end subroutine
      end interface

      interface
        subroutine dfft3f (n1, n2, n3, a, lda, mda, b, ldb,mdb)
          integer    n1, n2, n3, lda, mda, ldb, mdb
          complex    *16 a(*), b(*)
        end subroutine
      end interface

      interface
        subroutine df2t3f (n1, n2, n3, a, lda, mda, b, ldb,mdb, wff1, wff2, wff3, cpy)
          integer    lda, mda, ldb, mdb, n1, n2, n3
          double precision wff1(*), wff2(*), wff3(*), cpy(*)
          complex    *16 a(lda,mda,*), b(lda,mda,*)
        end subroutine
      end interface

      interface
        subroutine dfft3b (n1, n2, n3, a, lda, mda, b, ldb,mdb)
          integer    n1, n2, n3, lda, mda, ldb, mdb
          complex    *16 a(*), b(*)
        end subroutine
      end interface

      interface
        subroutine df2t3b (n1, n2, n3, a, lda, mda, b, ldb,mdb, wff1, wff2, wff3, cpy)
          integer    n1, n2, n3, lda, mda, ldb, mdb
          double precision wff1(*), wff2(*), wff3(*), cpy(*)
          complex    *16 a(lda,mda,*), b(lda,mda,*)
        end subroutine
      end interface

      interface
        subroutine drconv (ido, nx, x, ny, y, ipad, nz, z,zhat)
          integer    ido, nx, ny, ipad, nz
          double precision x(*), y(*), z(*), zhat(*)
        end subroutine
      end interface

      interface
        subroutine dr2onv (ido, nx, x, ny, y, ipad, nz, z,zhat, xwk, ywk, wk)
          integer    ido, nx, ny, ipad, nz
          double precision x(*), y(*), z(*), zhat(*), xwk(*), ywk(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dcconv (ido, nx, x, ny, y, ipad, nz, z,zhat)
          integer    ido, nx, ny, ipad, nz
          complex    *16 x(*), y(*), z(*), zhat(*)
        end subroutine
      end interface

      interface
        subroutine dc2onv (ido, nx, x, ny, y, ipad, nz, z,zhat, xwk, ywk, wk)
          integer    ido, nx, ny, ipad, nz
          double precision wk(*)
          complex    *16 x(*), y(*), z(*), zhat(*), xwk(*), ywk(*)
        end subroutine
      end interface

      interface
        subroutine drcorl (ido, n, x, y, ipad, nz, z, zhat)
          integer    ido, n, ipad, nz
          double precision x(*), y(*), z(*), zhat(*)
        end subroutine
      end interface

      interface
        subroutine dr2orl (ido, n, x, y, ipad, nz, z, zhat,xwk, ywk, wk)
          integer    n, ipad, ido, nz
          double precision x(*), y(*), z(*), zhat(*), xwk(*), ywk(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dccorl (ido, n, x, y, ipad, nz, z, zhat)
          integer    ido, n, ipad, nz
          complex    *16 x(*), y(*), z(*), zhat(*)
        end subroutine
      end interface

      interface
        subroutine dc2orl (ido, n, x, y, ipad, nz, z, zhat,xwk, ywk, wk)
          integer    ido, n, ipad, nz
          double precision wk(*)
          complex    *16 x(*), y(*), z(*), zhat(*), xwk(*), ywk(*)
        end subroutine
      end interface

      interface
        subroutine dinlap (f, n, t, alpha, relerr, kmax,finv)
          integer    n, kmax
          double precision alpha, relerr, t(*), finv(*)
          complex    *16 f
          external   f
        end subroutine
      end interface

      interface
        subroutine dsinlp (f, n, t, sigma0, epstol, errvec,finv)
          integer    n
          double precision sigma0, epstol, t(*), errvec(*), finv(*)
          complex    *16 f
          external   f
        end subroutine
      end interface

      interface
        subroutine ds2nlp (f, n, t, sigma0, epstol, errvec,finv, sigma, bvalue, mtop, wk, iflovc)
          integer    n, mtop, iflovc(*)
          double precision sigma0, epstol, sigma, bvalue, t(*), errvec(*),  &
     &           finv(*), wk(*)
          complex    *16 f
          external   f
        end subroutine
      end interface

!
!     Chapter 7:  Nonlinear Equations
!

      interface
        subroutine dzplrc (ndeg, coeff, root)
          integer    ndeg
          double precision coeff(*)
          complex    *16 root(*)
        end subroutine
      end interface

      interface
        subroutine dzporc (ndeg, coeff, root)
          integer    ndeg
          double precision coeff(*)
          complex    *16 root(*)
        end subroutine
      end interface

      interface
        subroutine dzpocc (ndeg, coeff, root)
          integer    ndeg
          complex    *16 coeff(*), root(*)
        end subroutine
      end interface

      interface
        subroutine dzanly (f, errabs, errrel, nknown, nnew,nguess, xinit, itmax, x, infer)
          integer    nknown, nnew, nguess, itmax, infer(*)
          double precision errabs, errrel
          complex    *16 f, xinit(*), x(*)
          external   f
        end subroutine
      end interface

      interface
        subroutine dzbren (f, errabs, errrel, a, b, maxfn)
          integer    maxfn
          double precision f, errabs, errrel, a, b
          external   f
        end subroutine
      end interface

      interface
        subroutine dzreal (f, errabs, errrel, eps, eta,nroot, itmax, xguess, x, infer)
          integer    nroot, itmax, infer(*)
          double precision f, errabs, errrel, eps, eta, xguess(*), x(*)
          external   f
        end subroutine
      end interface

      interface
        subroutine dneqnf (fcn, errrel, n, itmax, xinit, x,fnorm)
          integer    n, itmax
          double precision errrel, fnorm, xinit(*), x(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dn2qnf (fcn, errrel, n, itmax, xinit, x,fnorm, fvec, fjac, r, qtf, wk)
          integer    n, itmax
          double precision errrel, fnorm, xinit(*), x(*), fvec(*),          &
     &           fjac(n,*), r(*), qtf(*), wk(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dneqnj (fcn, lsjac, errrel, n, itmax,xinit, x, fnorm)
          integer    n, itmax
          double precision errrel, fnorm, xinit(*), x(*)
          external   fcn, lsjac
        end subroutine
      end interface

      interface
        subroutine dn2qnj (fcn, lsjac, errrel, n, itmax,xinit, x, fnorm, fvec, fjac, r, qtf, wk)
          integer    n, itmax
          double precision errrel, fnorm, xinit(*), x(*), fvec(*),          &
     &           fjac(n,*), r(*), qtf(*), wk(*)
          external   fcn, lsjac
        end subroutine
      end interface

      interface
        subroutine dneqbf (fcn, n, xguess, xscale, fscale,iparam, rparam, x, fvec)
          integer    n, iparam(*)
          double precision xguess(*), xscale(*), fscale(*), rparam(*),      &
     &           x(*), fvec(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dn2qbf (fcn, n, xguess, xscale, fscale,iparam, rparam, x, fvec, wk, lwk)
          integer    n, lwk, iparam(*)
          double precision xguess(*), xscale(*), fscale(*), rparam(*),      &
     &           x(*), fvec(*), wk(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dneqbj (fcn, jac, n, xguess, xscale,fscale, iparam, rparam, x, fvec)
          integer    n, iparam(*)
          double precision xguess(*), xscale(*), fscale(*), rparam(*),      &
     &           x(*), fvec(*)
          external   fcn, jac
        end subroutine
      end interface

      interface
        subroutine dn2qbj (fcn, jac, n, xguess, xscale,fscale, iparam, rparam, x, fvec, wk, lwk)
          integer    n, lwk, iparam(*)
          double precision xguess(*), xscale(*), fscale(*), rparam(*),      &
     &           x(*), fvec(*), wk(*)
          external   fcn, jac
        end subroutine
      end interface

      interface
        subroutine dn4qbj (iparam, rparam)
          integer    iparam(*)
          double precision rparam(*)
        end subroutine
      end interface

!
!     Chapter 8:  Optimization
!

      interface
        subroutine duvmif (f, xguess, step, bound, xacc,maxfn, x)
          integer    maxfn
          double precision f, xguess, step, bound, xacc, x
          external   f
        end subroutine
      end interface

      interface
        subroutine duvmid (f, g, xguess, errrel, gtol, maxfn,a, b, x, fx, gx)
          integer    maxfn
          double precision f, g, xguess, errrel, gtol, a, b, x, fx, gx
          external   f, g
        end subroutine
      end interface

      interface
        subroutine duvmgs (f, a, b, tol, xmin)
          double precision f, a, b, tol, xmin
          external   f
        end subroutine
      end interface

      interface
        subroutine duming (fcn, grad, n, xguess, xscale,fscale, iparam, rparam, x, fvalue)
          integer    n, iparam(*)
          double precision fscale, fvalue, xguess(*), xscale(*),            &
     &           rparam(*), x(*)
          external   fcn, grad
        end subroutine
      end interface

      interface
        subroutine du2ing (fcn, grad, n, xguess, xscale,fscale, iparam, rparam, x, fvalue, wk)
          integer    n, iparam(*)
          double precision xguess(*), fscale, fvalue, xscale(*),            &
     &           rparam(*), x(*), wk(*)
          external   fcn, grad
        end subroutine
      end interface

      interface
        subroutine dumidh (fcn, grad, n, xguess, xscale,fscale, iparam, rparam, x, fvalue)
          integer    n, iparam(*)
          double precision fscale, fvalue, xguess(*), xscale(*),            &
     &           rparam(*), x(*)
          external   fcn, grad
        end subroutine
      end interface

      interface
        subroutine du2idh (fcn, grad, n, xguess, xscale,fscale, iparam, rparam, x, fvalue, wk)
          integer    n, iparam(*)
          double precision fscale, fvalue, xguess(*), xscale(*),            &
     &           rparam(*), x(*), wk(*)
          external   fcn, grad
        end subroutine
      end interface

      interface
        subroutine dumiah (fcn, grad, hess, n, xguess,xscale, fscale, iparam, rparam, x, fvalue)
          integer    n, iparam(*)
          double precision fscale, fvalue, xguess(*), xscale(*),            &
     &           rparam(*), x(*)
          external   fcn, grad, hess
        end subroutine
      end interface

      interface
        subroutine du2iah (fcn, grad, hess, n, xguess,xscale, fscale, iparam, rparam, x, fvalue, wk)
          integer    n, iparam(*)
          double precision fscale, fvalue, xguess(*), xscale(*),            &
     &           rparam(*), x(*), wk(*)
          external   fcn, grad, hess
        end subroutine
      end interface

      interface
        subroutine dumcgf (fcn, n, xguess, xscale, gradtl,maxfn, dfpred, x, g, fvalue)
          integer    n, maxfn
          double precision gradtl, dfpred, fvalue, xguess(*), xscale(*),    &
     &           x(*), g(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine du2cgf (fcn, n, xguess, xscale, gradtl,maxfn, dfpred, x, g, fvalue, s, rss, rsg, ginit, xopt,gopt)
          integer    n, maxfn
          double precision gradtl, dfpred, fvalue, xguess(*), xscale(*),    &
     &           x(*), g(*), s(*), rss(*), rsg(*), ginit(*), xopt(*),   &
     &           gopt(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dumcgg (fcn, grad, n, xguess, gradtl,maxfn, dfpred, x, g, fvalue)
          integer    n, maxfn
          double precision gradtl, dfpred, fvalue, xguess(*), x(*), g(*)
          external   fcn, grad
        end subroutine
      end interface

      interface
        subroutine du2cgg (fcn, grad, n, xguess, gradtl,maxfn, dfpred, x, g, fvalue, s, rss, rsg, ginit, xopt,gopt)
          integer    n, maxfn
          double precision gradtl, dfpred, fvalue, xguess(*), x(*), g(*),   &
     &           s(*), rss(*), rsg(*), ginit(*), xopt(*), gopt(*)
          external   fcn, grad
        end subroutine
      end interface

      interface
        subroutine dumpol (fcn, n, xguess, s, ftol, maxfcn,x, fvalue)
          integer    n, maxfcn
          double precision s, ftol, fvalue, xguess(*), x(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine du2pol (fcn, n, xguess, s, ftol, maxfcn,x, fvalue, wk)
          integer    n, maxfcn
          double precision s, ftol, fvalue, xguess(*), x(*), wk(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dunlsf (fcn, m, n, xguess, xscale, fscale,iparam, rparam, x, fvec, fjac, ldfjac)
          integer    m, n, ldfjac, iparam(*)
          double precision xguess(*), xscale(*), fscale(*), rparam(*),      &
     &           x(*), fvec(*), fjac(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine du2lsf (fcn, m, n, xguess, xscale, fscale,iparam, rparam, x, fvec, fjac, ldfjac, wk, iwk)
          integer    m, n, ldfjac, iparam(*), iwk(*)
          double precision xguess(*), xscale(*), fscale(*), rparam(*),      &
     &           x(*), fvec(*), fjac(*), wk(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine du4lsf (iparam, rparam)
          integer    iparam(*)
          double precision rparam(*)
        end subroutine
      end interface

      interface
        subroutine dunlsj (fcn, jac, m, n, xguess, xscale,fscale, iparam, rparam, x, fvec, fjac, ldfjac)
          integer    m, n, ldfjac, iparam(*)
          double precision xguess(*), xscale(*), fscale(*), rparam(*),      &
     &           x(*), fvec(*), fjac(ldfjac,*)
          external   fcn, jac
        end subroutine
      end interface

      interface
        subroutine du2lsj (fcn, jac, m, n, xguess, xscale,fscale, iparam, rparam, x, fvec, fjac, ldfjac, wk, iwk)
          integer    m, n, ldfjac, iparam(*), iwk(*)
          double precision xguess(*), xscale(*), fscale(*), rparam(*),      &
     &           x(*), fvec(*), fjac(ldfjac,*), wk(*)
          external   fcn, jac
        end subroutine
      end interface

      interface
        subroutine dbconf (fcn, n, xguess, ibtype, xlb, xub,xscale, fscale, iparam, rparam, x, fvalue)
          integer    n, ibtype, iparam(*)
          double precision fscale, fvalue, xguess(*), xlb(*), xub(*),       &
     &           xscale(*), rparam(*), x(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine db2onf (fcn, n, xguess, ibtype, xlb, xub,xscale, fscale, iparam, rparam, x, fvalue, wk, iwk)
          integer    n, ibtype, iparam(*), iwk(*)
          double precision fscale, fvalue, xguess(*), xlb(*), xub(*),       &
     &           xscale(*), rparam(*), x(*), wk(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dbcong (fcn, grad, n, xguess, ibtype, xlb,xub, xscale, fscale, iparam, rparam, x, fvalue)
          integer    n, ibtype, iparam(*)
          double precision fscale, fvalue, xguess(*), xlb(*), xub(*),       &
     &           xscale(*), rparam(*), x(*)
          external   fcn, grad
        end subroutine
      end interface

      interface
        subroutine db2ong (fcn, grad, n, xguess, ibtype, xlb,xub, xscale, fscale, iparam, rparam, x, fvalue, wk,iwk)
          integer    n, ibtype, iparam(*), iwk(*)
          double precision fscale, fvalue, xguess(*), xlb(*), xub(*),       &
     &           xscale(*), rparam(*), x(*), wk(*)
          external   fcn, grad
        end subroutine
      end interface

      interface
        subroutine dbcodh (fcn, grad, n, xguess, ibtype, xlb,xub, xscale, fscale, iparam, rparam, x, fvalue)
          integer    n, ibtype, iparam(*)
          double precision fscale, fvalue, xguess(*), xlb(*), xub(*),       &
     &           xscale(*), rparam(*), x(*)
          external   fcn, grad
        end subroutine
      end interface

      interface
        subroutine db2odh (fcn, grad, n, xguess, ibtype, xlb,xub, xscale, fscale, iparam, rparam, x, fvalue, wk,iwk)
          integer    n, ibtype, iparam(*), iwk(*)
          double precision fscale, fvalue, xguess(*), xlb(*), xub(*),       &
     &           xscale(*), rparam(*), x(*), wk(*)
          external   fcn, grad
        end subroutine
      end interface

      interface
        subroutine dbcoah (fcn, grad, hess, n, xguess,ibtype, xlb, xub, xscale, fscale, iparam, rparam, x,fvalue)
          integer    n, ibtype, iparam(*)
          double precision fscale, fvalue, xguess(*), xlb(*), xub(*),       &
     &           xscale(*), rparam(*), x(*)
          external   fcn, grad, hess
        end subroutine
      end interface

      interface
        subroutine db2oah (fcn, grad, hess, n, xguess,ibtype, xlb, xub, xscale, fscale, iparam, rparam, x,fvalue, wk, iwk)
          integer    n, ibtype, iparam(*), iwk(*)
          double precision fscale, fvalue, xguess(*), xlb(*), xub(*),       &
     &           xscale(*), rparam(*), x(*), wk(*)
          external   fcn, grad, hess
        end subroutine
      end interface

      interface
        subroutine dbcpol (fcn, n, xguess, ibtype, xlb, xub,ftol, maxfcn, x, fvalue)
          integer    n, ibtype, maxfcn
          double precision ftol, fvalue, xguess(*), xlb(*), xub(*), x(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine db2pol (fcn, n, xguess, ibtype, xlb, xub,ftol, maxfcn, x, fvalue, wk)
          integer    n, ibtype, maxfcn
          double precision ftol, fvalue, xguess(*), xlb(*), xub(*), x(*),   &
     &           wk(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dbclsf (fcn, m, n, xguess, ibtype, xlb,xub, xscale, fscale, iparam, rparam, x, fvec, fjac,ldfjac)
          integer    m, n, ibtype, ldfjac, iparam(*)
          double precision xguess(*), xlb(*), xub(*), xscale(*),            &
     &           fscale(*), rparam(*), x(*), fvec(*), fjac(ldfjac,*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine db2lsf (fcn, m, n, xguess, ibtype, xlb,xub, xscale, fscale, iparam, rparam, x, fvec, fjac,ldfjac, wk, iwk)
          integer    m, n, ibtype, ldfjac, iparam(*), iwk(*)
          double precision xguess(*), xlb(*), xub(*), xscale(*),            &
     &           fscale(*), rparam(*), x(*), fvec(*), fjac(*), wk(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dbclsj (fcn, jac, m, n, xguess, ibtype,xlb, xub, xscale, fscale, iparam, rparam, x, fvec,fjac, ldfjac)
          integer    m, n, ibtype, ldfjac, iparam(*)
          double precision xguess(*), xlb(*), xub(*), xscale(*),            &
     &           fscale(*), rparam(*), x(*), fvec(*), fjac(ldfjac,*)
          external   fcn, jac
        end subroutine
      end interface

      interface
        subroutine db2lsj (fcn, jac, m, n, xguess, ibtype,xlb, xub, xscale, fscale, iparam, rparam, x, fvec,fjac, ldfjac, wk, iwk)
          integer    m, n, ibtype, ldfjac, iparam(*), iwk(*)
          double precision xguess(*), xlb(*), xub(*), xscale(*),            &
     &           fscale(*), rparam(*), x(*), fvec(*), fjac(ldfjac,*),   &
     &           wk(*)
          external   fcn, jac
        end subroutine
      end interface

      interface
        subroutine ddlprs (m, nvar, a, lda, bl, bu, c,irtype, xlb, xub, obj, xsol, dsol)
          integer    m, nvar, lda, irtype(*)
          double precision obj, a(lda,*), bl(*), bu(*), c(*), xlb(*),       &
     &           xub(*), xsol(*), dsol(*)
        end subroutine
      end interface

      interface
        subroutine dd2prs (m, nvar, a, lda, bl, bu, c,irtype, xlb, xub, obj, xsol, dsol, awk, ldawk, wk, iwk)
          integer    m, nvar, lda, ldawk, irtype(*), iwk(*)
          double precision obj, a(lda,*), bl(*), bu(*), c(*), xlb(*),       &
     &           xub(*), xsol(*), dsol(*), awk(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dqprog (nvar, ncon, neq, a, lda, b, g, h,ldh, diag, sol, nact, iact, alamda)
          integer    nvar, ncon, neq, lda, ldh, nact, iact(*)
          double precision diag, a(lda,*), b(*), g(*), h(ldh,*), sol(*),    &
     &           alamda(*)
        end subroutine
      end interface

      interface
        subroutine dq2rog (nvar, ncon, neq, a, lda, b, grad,h, ldh, diag, sol, nact, iact, alamda, wk)
          integer    nvar, ncon, neq, lda, ldh, nact, iact(*)
          double precision diag, a(lda,*), b(*), grad(*), h(ldh,*),         &
     &           sol(*), alamda(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dlconf (fcn, nvar, ncon, neq, a, lda, b,xlb, xub, xguess, acc, maxfcn, sol, obj, nact, iact,alamda)
          integer    nvar, ncon, neq, lda, maxfcn, nact, iact(*)
          double precision acc, obj, a(lda,*), b(*), xlb(*), xub(*),        &
     &           xguess(*), sol(*), alamda(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dl2onf (fcn, n, m, meq, a, ia, b, xl, xu,x, acc, maxfcn, sol, obj, nact, iact, par, iprint,info, w)
          integer    n, m, meq, ia, maxfcn, nact, iprint, info, iact(*)
          double precision acc, obj, a(ia,*), b(*), xl(*), xu(*), x(*),     &
     &           sol(*), par(*), w(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dlcong (fcn, grad, nvar, ncon, neq, a,lda, b, xlb, xub, xguess, acc, maxfcn, sol, obj, nact,iact, alamda)
          integer    nvar, ncon, neq, lda, maxfcn, nact, iact(*)
          double precision acc, obj, a(lda,*), b(*), xlb(*), xub(*),        &
     &           xguess(*), sol(*), alamda(*)
          external   fcn, grad
        end subroutine
      end interface

      interface
        subroutine dl2ong (fcn, grad, n, m, meq, a, ia, b,xl, xu, x, acc, maxfcn, sol, obj, nact, iact, par,iprint, info, w)
          integer    n, m, meq, ia, maxfcn, nact, iprint, info, iact(*)
          double precision acc, obj, a(ia,*), b(*), xl(*), xu(*), x(*),     &
     &           sol(*), par(*), w(*)
          external   fcn, grad
        end subroutine
      end interface

      interface
        subroutine dnconf (fcn, m, me, n, xguess, ibtype,xlb, xub, xscale, iprint, maxitn, x, fvalue)
          integer    m, me, n, ibtype, iprint, maxitn
          double precision fvalue, xguess(*), xlb(*), xub(*), xscale(*),    &
     &           x(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dn2onf (fcns, m, me, n, xguess, ibtype,xlb, xub, xscale, iprint, maxitn, x, fvalue, wk, lwk,iwk, liwk, conwk)
          integer    m, me, n, ibtype, iprint, maxitn, lwk, liwk, iwk(*)
          double precision fvalue, xguess(*), xlb(*), xub(*), xscale(*),    &
     &           x(*), wk(*), conwk(*)
          external   fcns
        end subroutine
      end interface

      interface
        subroutine dn0onf (ido, m, me, n, ibtype, xlb, xub,iprint, maxitn, x, fvalue, g, df, dg, lddg, u, c, ldc,d, acc, scbou, maxfun, active, mode, wk, iwk, conwk)
          integer    ido, m, me, n, ibtype, iprint, maxitn, lddg, ldc,      &
     &           maxfun, mode, iwk(*)
          double precision fvalue, acc, scbou, xlb(*), xub(*), x(*), g(*),  &
     &           df(*), dg(lddg,*), u(*), c(ldc,*), d(*), wk(*),        &
     &           conwk(*)
          logical    active(*)
        end subroutine
      end interface

      interface
        subroutine dncong (fcn, grad, m, me, n, xguess,ibtype, xlb, xub, iprint, maxitn, x, fvalue)
          integer    m, me, n, ibtype, iprint, maxitn
          double precision fvalue, xguess(*), xlb(*), xub(*), x(*)
          external   fcn, grad
        end subroutine
      end interface

      interface
        subroutine dn2ong (fcns, grad, m, me, n, xguess,ibtype, xlb, xub, iprint, maxitn, x, fvalue, wk, lwk,iwk, liwk)
          integer    m, me, n, ibtype, iprint, maxitn, lwk, liwk, iwk(*)
          double precision fvalue, xguess(*), xlb(*), xub(*), x(*), wk(*)
          external   fcns, grad
        end subroutine
      end interface

      interface
        subroutine dcdgrd (fcn, n, xc, xscale, epsfcn, gc)
          integer    n
          double precision epsfcn, xc(*), xscale(*), gc(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dfdgrd (fcn, n, xc, xscale, fc, epsfcn,gc)
          integer    n
          double precision fc, epsfcn, xc(*), xscale(*), gc(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dfdhes (fcn, n, xc, xscale, fc, epsfcn, h,ldh)
          integer    n, ldh
          double precision fc, epsfcn, xc(*), xscale(*), h(ldh,*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine df2hes (fcn, n, xc, xscale, fc, epsfcn, h,ldh, stepsz, fneibr)
          integer    n, ldh
          double precision fc, epsfcn, xc(*), xscale(*), h(ldh,*),          &
     &           stepsz(*), fneibr(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dgdhes (grad, n, xc, xscale, gc, epsfcn,h, ldh)
          integer    n, ldh
          double precision epsfcn, xc(*), xscale(*), gc(*), h(ldh,*)
          external   grad
        end subroutine
      end interface

      interface
        subroutine dg2hes (grad, n, xc, xscale, gc, epsfcn,h, ldh, wk)
          integer    n, ldh
          double precision epsfcn, xc(*), xscale(*), gc(*), h(ldh,*), wk(*)
          external   grad
        end subroutine
      end interface

      interface
        subroutine dfdjac (fcn, m, n, xc, xscale, fc, epsfcn,fjac, ldfjac)
          integer    m, n, ldfjac
          double precision epsfcn, xc(*), xscale(*), fc(*), fjac(ldfjac,*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine df2jac (fcn, m, n, xc, xscale, fc, epsfcn,fjac, ldfjac, work)
          integer    m, n, ldfjac
          double precision epsfcn, xc(*), xscale(*), fc(*),                 &
     &           fjac(ldfjac,*), work(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dchgrd (fcn, grad, n, x, info)
          integer    n, info(*)
          double precision grad(*), x(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dc2grd (fcn, grad, n, x, info, fx, xscale,epsfcn, xnew)
          integer    n, info(*)
          double precision fx, epsfcn, grad(*), x(*), xscale(*), xnew(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine dchhes (grad, hess, n, x, info, ldinfo)
          integer    n, ldinfo, info(ldinfo,*)
          double precision x(*)
          external   grad, hess
        end subroutine
      end interface

      interface
        subroutine dc2hes (grad, hess, n, x, info, ldinfo,gx, hx, hs, xscale, epsfcn, inft, xnew)
          integer    n, ldinfo, info(ldinfo,*), inft(*)
          double precision epsfcn, x(*), gx(*), hx(n,*), hs(*), xscale(*),  &
     &           xnew(*)
          external   grad, hess
        end subroutine
      end interface

      interface
        subroutine dchjac (fcn, jac, m, n, x, info, ldinfo)
          integer    m, n, ldinfo, info(ldinfo,*)
          double precision x(*)
          external   fcn, jac
        end subroutine
      end interface

      interface
        subroutine dc2jac (fcn, jac, m, n, x, info, ldinfo,fx, fjac, grad, xscale, epsfcn, inft, xnew)
          integer    m, n, ldinfo, info(ldinfo,*), inft(*)
          double precision epsfcn, x(*), fx(*), fjac(m,*), grad(*),         &
     &           xscale(*), xnew(*)
          external   fcn, jac
        end subroutine
      end interface

      interface
        subroutine dggues (n, a, b, k, ido, s)
          integer    n, k, ido
          double precision a(*), b(*), s(*)
        end subroutine
      end interface

      interface
        subroutine dg2ues (n, a, b, k, ido, s, wk, iwk)
          integer    n, k, ido, iwk(*)
          double precision a(*), b(*), s(*), wk(*)
        end subroutine
      end interface

!
!     Chapter 9:  Basic Matrix/Vector Operations
!

      interface
        subroutine dcrgrg (n, a, lda, b, ldb)
          integer    n, lda, ldb
          double precision a(lda,*), b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine dccgcg (n, a, lda, b, ldb)
          integer    n, lda, ldb
          complex    *16 a(lda,*), b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine dcrbrb (n, a, lda, nlca, nuca, b, ldb,nlcb, nucb)
          integer    n, lda, nlca, nuca, ldb, nlcb, nucb
          double precision a(lda,*), b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine dccbcb (n, a, lda, nlca, nuca, b, ldb,nlcb, nucb)
          integer    n, lda, nlca, nuca, ldb, nlcb, nucb
          complex    *16 a(lda,*), b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine dcrgrb (n, a, lda, nlc, nuc, b, ldb)
          integer    n, lda, nlc, nuc, ldb
          double precision a(lda,*), b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine dcrbrg (n, a, lda, nlc, nuc, b, ldb)
          integer    n, lda, nlc, nuc, ldb
          double precision a(lda,*), b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine dccgcb (n, a, lda, nlc, nuc, b, ldb)
          integer    n, lda, nlc, nuc, ldb
          complex    *16 a(lda,*), b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine dccbcg (n, a, lda, nlc, nuc, b, ldb)
          integer    n, lda, nlc, nuc, ldb
          complex    *16 a(lda,*), b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine dcrgcg (n, a, lda, b, ldb)
          integer    n, lda, ldb
          double precision a(lda,*)
          complex    *16 b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine dcrrcr (nra, nca, a, lda, nrb, ncb, b,ldb)
          integer    nra, nca, lda, nrb, ncb, ldb
          double precision a(lda,*)
          complex    *16 b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine dcrbcb (n, a, lda, nlca, nuca, b, ldb,nlcb, nucb)
          integer    n, lda, nlca, nuca, ldb, nlcb, nucb
          double precision a(lda,*)
          complex    *16 b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine dcsfrg (n, a, lda)
          integer    n, lda
          double precision a(lda,*)
        end subroutine
      end interface

      interface
        subroutine dchfcg (n, a, lda)
          integer    n, lda
          complex    *16 a(lda,*)
        end subroutine
      end interface

      interface
        subroutine dcsbrb (n, a, lda, nuca, b, ldb, nlcb,nucb)
          integer    n, lda, nuca, ldb, nlcb, nucb
          double precision a(lda,*), b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine dchbcb (n, a, lda, nuca, b, ldb, nlcb,nucb)
          integer    n, lda, nuca, ldb, nlcb, nucb
          complex    *16 a(lda,*), b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine dmxytf (nra, nca, a, lda, nrb, ncb, b,ldb, nrc, ncc, c, ldc)
          integer    nra, nca, lda, nrb, ncb, ldb, nrc, ncc, ldc
          double precision a(lda,*), b(ldb,*), c(ldc,*)
        end subroutine
      end interface

      interface
        subroutine dmcrcr (nra, nca, a, lda, nrb, ncb, b,ldb, nrc, ncc, c, ldc)
          integer    nra, nca, lda, nrb, ncb, ldb, nrc, ncc, ldc
          complex    *16 a(lda,*), b(ldb,*), c(ldc,*)
        end subroutine
      end interface

      interface
        subroutine dhrrrr (nra, nca, a, lda, nrb, ncb, b,ldb, nrc, ncc, c, ldc)
          integer    nra, nca, lda, nrb, ncb, ldb, nrc, ncc, ldc
          double precision a(lda,*), b(ldb,*), c(ldc,*)
        end subroutine
      end interface

      interface
        subroutine dpolrg (n, a, lda, ncoef, coef, b, ldb)
          integer    n, lda, ncoef, ldb
          double precision a(lda,*), coef(*), b(ldb,*)
        end subroutine
      end interface

      interface
        subroutine dp2lrg (n, a, lda, ncoef, coef, b, ldb,work)
          integer    n, lda, ncoef, ldb
          double precision a(lda,*), coef(*), b(ldb,*), work(n,*)
        end subroutine
      end interface

      interface
        subroutine dmurrv (nra, nca, a, lda, nx, x, ipath,ny, y)
          integer    nra, nca, lda, nx, ipath, ny
          double precision a(lda,*), x(*), y(*)
        end subroutine
      end interface

      interface
        subroutine dmurbv (n, a, lda, nlca, nuca, nx, x,ipath, ny, y)
          integer    n, lda, nlca, nuca, nx, ipath, ny
          double precision a(lda,n), x(nx), y(ny)
        end subroutine
      end interface

      interface
        subroutine dmucrv (nra, nca, a, lda, nx, x, ipath,ny, y)
          integer    nra, nca, lda, nx, ipath, ny
          complex    *16 a(lda,*), x(*), y(*)
        end subroutine
      end interface

      interface
        subroutine dmucbv (n, a, lda, nlca, nuca, nx, x,ipath, ny, y)
          integer    n, lda, nlca, nuca, nx, ipath, ny
          complex    *16 a(lda,n), x(nx), y(ny)
        end subroutine
      end interface

      interface
        subroutine darbrb (n, a, lda, nlca, nuca, b, ldb,nlcb, nucb, c, ldc, nlcc, nucc)
          integer    n, lda, nlca, nuca, ldb, nlcb, nucb, ldc, nlcc, nucc
          double precision a(lda,*), b(ldb,*), c(ldc,*)
        end subroutine
      end interface

      interface
        subroutine dacbcb (n, a, lda, nlca, nuca, b, ldb,nlcb, nucb, c, ldc, nlcc, nucc)
          integer    n, lda, nlca, nuca, ldb, nlcb, nucb, ldc, nlcc, nucc
          complex    *16 a(lda,*), b(ldb,*), c(ldc,*)
        end subroutine
      end interface

      interface
        subroutine dnrirr (nra, nca, a, lda, anorm)
          integer    nra, nca, lda
          double precision anorm, a(lda,*)
        end subroutine
      end interface

      interface
        subroutine dnr1rr (nra, nca, a, lda, anorm)
          integer    nra, nca, lda
          double precision anorm, a(lda,*)
        end subroutine
      end interface

      interface
        subroutine dnr2rr (nra, nca, a, lda, anorm)
          integer    nra, nca, lda
          double precision anorm, a(lda,*)
        end subroutine
      end interface

      interface
        subroutine dnr1rb (n, a, lda, nlca, nuca, anorm)
          integer    n, lda, nlca, nuca
          double precision anorm, a(lda,*)
        end subroutine
      end interface

      interface
        subroutine dnr1cb (n, a, lda, nlca, nuca, anorm)
          integer    n, lda, nlca, nuca
          double precision anorm
          complex    *16 a(lda,*)
        end subroutine
      end interface

      interface
        double precision function ddisl2 (n, x, incx, y,incy)
          integer    n, incx, incy
          double precision x(*), y(*)
        end function
      end interface

      interface
        double precision function ddisl1 (n, x, incx, y,incy)
          integer    n, incx, incy
          double precision x(*), y(*)
        end function
      end interface

      interface
        double precision function ddisli (n, x, incx, y,incy)
          integer    n, incx, incy
          double precision x(*), y(*)
        end function
      end interface

      interface
        subroutine dvconr (nx, x, ny, y, nz, z)
          integer    nx, ny, nz
          double precision x(*), y(*), z(*)
        end subroutine
      end interface

      interface
        subroutine dv2onr (nx, x, ny, y, nz, z, xwk, ywk,zwk, wk)
          integer    nx, ny, nz
          double precision x(*), y(*), z(*), wk(*)
          complex    *16 xwk(*), ywk(*), zwk(*)
        end subroutine
      end interface

      interface
        subroutine dvconc (nx, x, ny, y, nz, z)
          integer    nx, ny, nz
          complex    *16 x(*), y(*), z(*)
        end subroutine
      end interface

      interface
        subroutine dv2onc (nx, x, ny, y, nz, z, xwk, ywk, wk)
          integer    nx, ny, nz
          double precision wk(*)
          complex    *16 x(*), y(*), z(*), xwk(*), ywk(*)
        end subroutine
      end interface

!
!     Chapter 10:  Utilities
!

      interface
        subroutine dwrcrn (title, nra, nca, a, lda, itring)
          integer    nra, nca, lda, itring
          complex    *16 a(*)
          character  title*(*)
        end subroutine
      end interface

      interface
        subroutine dwrcrl (title, nra, nca, a, lda, itring,fmt, rlabel, clabel)
          integer    nra, nca, lda, itring
          complex    *16 a(*)
          character  title*(*), fmt*(*), rlabel(*)*(*), clabel(0:*)*(*)
        end subroutine
      end interface

      interface
        subroutine dw2crl (title, nra, nca, a, lda, itring,fmt, rlabel, clabel, chwk)
          integer    nra, nca, lda, itring
          complex    *16 a(*)
          character  title*(*), fmt*(*), rlabel(*)*(*), clabel(0:*)*(*),    &
     &           chwk(*)*10
        end subroutine
      end interface

      interface
        subroutine dsvrbn (n, ra, rb)
          integer    n
          double precision ra(*), rb(*)
        end subroutine
      end interface

      interface
        subroutine dsvrbp (n, ra, rb, iperm)
          integer    n, iperm(*)
          double precision ra(*), rb(*)
        end subroutine
      end interface

      interface
        double precision function dconst (name)
          character  name*(*)
        end function
      end interface

      interface
        subroutine dcunit (x, xunits, y, yunits)
          double precision x, y
          character  xunits*(*), yunits*(*)
        end subroutine
      end interface

      interface
        double precision function dhypot (a, b)
          double precision a, b
        end function
      end interface

!
!     Special Functions Chapter 1:  Elementary Functions
!

      interface
        double precision function zarg (z)
          complex    *16 z
        end function
      end interface

      interface
        double precision function dcbrt (x)
          double precision x
        end function
      end interface

      interface
        complex *16 function zcbrt (z)
          complex    *16 z
        end function
      end interface

      interface
        double precision function dexprl (x)
          double precision x
        end function
      end interface

      interface
        complex *16 function zlog10 (z)
          complex    *16 z
        end function
      end interface

      interface
        double precision function dlnrel (x)
          double precision x
        end function
      end interface

!
!     Special Functions Chapter 2:  Trigonometric and Hyperbolic
!                                   Functions
!

      interface
        complex *16 function ztan (z)
          complex    *16 z
        end function
      end interface

      interface
        double precision function dcot (x)
          double precision x
        end function
      end interface

      interface
        complex *16 function zcot (z)
          complex    *16 z
        end function
      end interface

      interface
        double precision function dsindg (x)
          double precision x
        end function
      end interface

      interface
        double precision function dcosdg (x)
          double precision x
        end function
      end interface

      interface
        complex *16 function zasin (zinp)
          complex    *16 zinp
        end function
      end interface

      interface
        complex *16 function zacos (z)
          complex    *16 z
        end function
      end interface

      interface
        complex *16 function zatan (z)
          complex    *16 z
        end function
      end interface

      interface
        complex *16 function zatan2 (csn, ccs)
          complex    *16 csn, ccs
        end function
      end interface

      interface
        complex *16 function zsinh (z)
          complex    *16 z
        end function
      end interface

      interface
        complex *16 function zcosh (z)
          complex    *16 z
        end function
      end interface

      interface
        complex *16 function ztanh (z)
          complex    *16 z
        end function
      end interface

      interface
        double precision function dasinh (x)
          double precision x
        end function
      end interface

      interface
        complex *16 function zasinh (z)
          complex    *16 z
        end function
      end interface

      interface
        double precision function dacosh (x)
          double precision x
        end function
      end interface

      interface
        complex *16 function zacosh (z)
          complex    *16 z
        end function
      end interface

      interface
        double precision function datanh (x)
          double precision x
        end function
      end interface

      interface
        complex *16 function zatanh (z)
          complex    *16 z
        end function
      end interface

!
!     Special Functions Chapter 3:  Exponential Integrals and Related
!                                   Functions
!

      interface
        double precision function dei (x)
          double precision x
        end function
      end interface

      interface
        double precision function de1 (x)
          double precision x
        end function
      end interface

      interface
        subroutine dene (x, n, f)
          integer    n
          double precision x, f(*)
        end subroutine
      end interface

      interface
        double precision function dli (x)
          double precision x
        end function
      end interface

      interface
        double precision function dsi (x)
          double precision x
        end function
      end interface

      interface
        double precision function dci (x)
          double precision x
        end function
      end interface

      interface
        double precision function dcin (x)
          double precision x
        end function
      end interface

      interface
        double precision function dshi (x)
          double precision x
        end function
      end interface

      interface
        double precision function dchi (x)
          double precision x
        end function
      end interface

      interface
        double precision function dcinh (x)
          double precision x
        end function
      end interface

!
!     Special Functions Chapter 4:  Gamma Function and Related Functions
!

      interface
        double precision function dfac (n)
          integer    n
        end function
      end interface

      interface
        double precision function dgamr (x)
          double precision x
        end function
      end interface

      interface
        double precision function dlngam (x)
          double precision x
        end function
      end interface

      interface
        subroutine dlgams (x, algm, s)
          double precision x, algm, s
        end subroutine
      end interface

      interface
        double precision function dgami (a, x)
          double precision a, x
        end function
      end interface

      interface
        double precision function dgamic (a, x)
          double precision a, x
        end function
      end interface

      interface
        double precision function dgamit (a, x)
          double precision a, x
        end function
      end interface

      interface
        double precision function dpsi (x)
          double precision x
        end function
      end interface

      interface
        double precision function dpoch (a, x)
          double precision a, x
        end function
      end interface

      interface
        double precision function dpoch1 (a, x)
          double precision a, x
        end function
      end interface

      interface
        double precision function dlbeta (a, b)
          double precision a, b
        end function
      end interface

      interface
        double precision function dbetai (x, pin, qin)
          double precision x, pin, qin
        end function
      end interface

!
!     Special Functions Chapter 5:  Error Function and Related Functions
!

      interface
        double precision function derf (x)
          double precision x
        end function
      end interface

      interface
        double precision function derfc (x)
          double precision x
        end function
      end interface

      interface
        double precision function derfce (x)
          double precision x
        end function
      end interface

      interface
        complex *16 function zerfe (z)
          complex    *16 z
        end function
      end interface

      interface
        double precision function derfi (x)
          double precision x
        end function
      end interface

      interface
        double precision function derfci (x)
          double precision x
        end function
      end interface

      interface
        double precision function ddaws (x)
          double precision x
        end function
      end interface

      interface
        double precision function dfresc (x)
          double precision x
        end function
      end interface

      interface
        double precision function dfress (x)
          double precision x
        end function
      end interface

!
!     Special Functions Chapter 6:  Bessel Functions
!

      interface
        double precision function dbsj0 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbsj1 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbsy0 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbsy1 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbsi0 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbsi1 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbsk0 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbsk1 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbsi0e (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbsi1e (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbsk0e (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbsk1e (x)
          double precision x
        end function
      end interface

      interface
        subroutine dbsjns (x, n, bs)
          integer    n
          double precision x, bs(*)
        end subroutine
      end interface

      interface
        subroutine dcbjns (z, n, cbs)
          integer    n
          complex    *16 z, cbs(*)
        end subroutine
      end interface

      interface
        subroutine dbsins (x, n, bsi)
          integer    n
          double precision x, bsi(*)
        end subroutine
      end interface

      interface
        subroutine dcbins (z, n, cbs)
          integer    n
          complex    *16 z, cbs(*)
        end subroutine
      end interface

      interface
        subroutine dbsjs (xnu, x, n, bs)
          integer    n
          double precision xnu, x, bs(*)
        end subroutine
      end interface

      interface
        subroutine db2js (xnu, x, n, bs, wk)
          integer    n
          double precision xnu, x, bs(*), wk(*)
        end subroutine
      end interface

      interface
        subroutine dbsys (xnu, x, n, bsy)
          integer    n
          double precision xnu, x, bsy(*)
        end subroutine
      end interface

      interface
        subroutine dbsis (xnu, x, n, bsi)
          integer    n
          double precision xnu, x, bsi(*)
        end subroutine
      end interface

      interface
        subroutine dbsies (xnu, x, n, bsi)
          integer    n
          double precision xnu, x, bsi(*)
        end subroutine
      end interface

      interface
        subroutine dbsks (xnu, x, nin, bk)
          integer    nin
          double precision xnu, x, bk(*)
        end subroutine
      end interface

      interface
        subroutine dbskes (xnu, x, nin, bke)
          integer    nin
          double precision xnu, x, bke(*)
        end subroutine
      end interface

      interface
        subroutine dcbjs (xnu, z, n, cbs)
          integer    n
          double precision xnu
          complex    *16 z, cbs(*)
        end subroutine
      end interface

      interface
        subroutine dcbys (xnu, z, n, cbs)
          integer    n
          double precision xnu
          complex    *16 z, cbs(*)
        end subroutine
      end interface

      interface
        subroutine dc2ys (xnu, z, n, cbs, fk)
          integer    n
          double precision xnu
          complex    *16 z, cbs(*), fk(*)
        end subroutine
      end interface

      interface
        subroutine dcbis (xnu, z, n, cbs)
          integer    n
          double precision xnu
          complex    *16 z, cbs(*)
        end subroutine
      end interface

      interface
        subroutine dcbks (xnu, z, n, cbs)
          integer    n
          double precision xnu
          complex    *16 z, cbs(*)
        end subroutine
      end interface

      interface
        subroutine dc2ks (xnu, z, n, cbs, fi)
          integer    n
          double precision xnu
          complex    *16 z, cbs(*), fi(*)
        end subroutine
      end interface

!
!     Special Functions Chapter 7:  Kelvin Functions
!

      interface
        double precision function dber0 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbei0 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dker0 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dkei0 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dberp0 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbeip0 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dkerp0 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dkeip0 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dber1 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbei1 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dker1 (x)
          double precision x
        end function
      end interface

      interface
        double precision function dkei1 (x)
          double precision x
        end function
      end interface

!
!     Special Functions Chapter 8:  Airy Functions
!

      interface
        double precision function dai (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbi (x)
          double precision x
        end function
      end interface

      interface
        double precision function daid (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbid (x)
          double precision x
        end function
      end interface

      interface
        double precision function daie (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbie (x)
          double precision x
        end function
      end interface

      interface
        double precision function daide (x)
          double precision x
        end function
      end interface

      interface
        double precision function dbide (x)
          double precision x
        end function
      end interface

!
!     Special Functions Chapter 9:  Elliptic Integrals
!

      interface
        double precision function delk (x)
          double precision x
        end function
      end interface

      interface
        double precision function dele (x)
          double precision x
        end function
      end interface

      interface
        double precision function delrf (x, y, z)
          double precision x, y, z
        end function
      end interface

      interface
        double precision function delrd (x, y, z)
          double precision x, y, z
        end function
      end interface

      interface
        double precision function delrj (x, y, z, rho)
          double precision x, y, z, rho
        end function
      end interface

      interface
        double precision function delrc (x, y)
          double precision x, y
        end function
      end interface

!
!     Special Functions Chapter 10:  Elliptic and Related Functions
!

      interface
        complex *16 function zwpl (z)
          complex    *16 z
        end function
      end interface

      interface
        complex *16 function zwpld (z)
          complex    *16 z
        end function
      end interface

      interface
        complex *16 function zwpq (z)
          complex    *16 z
        end function
      end interface

      interface
        complex *16 function zwpqd (z)
          complex    *16 z
        end function
      end interface

      interface
        double precision function dejsn (x, am)
          double precision x, am
        end function
      end interface

      interface
        complex *16 function zejsn (z, am)
          double precision am
          complex    *16 z
        end function
      end interface

      interface
        double precision function dejcn (x, am)
          double precision x, am
        end function
      end interface

      interface
        complex *16 function zejcn (z, am)
          double precision am
          complex    *16 z
        end function
      end interface

      interface
        double precision function dejdn (x, am)
          double precision x, am
        end function
      end interface

      interface
        complex *16 function zejdn (z, am)
          double precision am
          complex    *16 z
        end function
      end interface

!
!     Special Functions Chapter 12:  Mathieu Functions
!

      interface
        subroutine dmatee (q, n, isym, iper, eval)
          integer    n, isym, iper
          double precision q, eval(*)
        end subroutine
      end interface

      interface
        subroutine dm2tee (q, n, isym, iper, eval, norder,workd, worke)
          integer    n, isym, iper, norder
          double precision q, eval(*), workd(*), worke(*)
        end subroutine
      end interface

      interface
        subroutine dm3tee (q, n, mmax)
          integer    n, mmax
          double precision q
        end subroutine
      end interface

      interface
        subroutine dmatce (x, q, n, ce)
          integer    n
          double precision x, q, ce(*)
        end subroutine
      end interface

      interface
        subroutine dm2tce (x, q, n, ce, norder, needev,eval0, eval1, coef, work, bsj)
          logical    needev
          integer    n, norder
          double precision x, q, ce(*), eval0(*), eval1(*), coef(*),        &
     &           work(*), bsj(*)
        end subroutine
      end interface

      interface
        subroutine dmatse (x, q, n, se)
          integer    n
          double precision x, q, se(*)
        end subroutine
      end interface

      interface
        subroutine dm2tse (x, q, n, se, norder, needev,eval0, eval1, coef, work, bsi)
          logical    needev
          integer    n, norder
          double precision x, q, se(*), eval0(*), eval1(*), coef(*),        &
     &           work(*), bsi(*)
        end subroutine
      end interface

!
!     Special Functions Chapter 13:  Miscellaneous Functions
!

      interface
        double precision function dspenc (x)
          double precision x
        end function
      end interface

      interface
        integer function initds (dos, nos, eta)
          integer    nos
          real       eta
          double precision dos(*)
        end function
      end interface

      interface
        double precision function dcsevl (x, cs, n)
          integer    n
          double precision x, cs(*)
        end function
      end interface

!
!     Deprecated Routines
!

      interface
        subroutine de2ahf (n, neval, a, lda, small, eval,acopy, rwk, cwk, iwk)
          integer    n, neval, lda, iwk(*)
          double precision eval(*), rwk(n,*)
          complex    *16 a(lda,*), acopy(n,*), cwk(n,*)
          logical    small
        end subroutine
      end interface

      interface
        subroutine de2asf (n, neval, a, lda, small, eval,acopy, wk, iwk)
          integer    n, neval, lda, iwk(*)
          double precision a(lda,*), eval(*), acopy(n,*), wk(*)
          logical    small
        end subroutine
      end interface

      interface
        subroutine de2bhf (n, mxeval, a, lda, elow, ehigh,neval, eval, acopy, rwk, cwk, iwk)
          integer    n, mxeval, lda, neval, iwk(*)
          double precision elow, ehigh, eval(*), rwk(n,*)
          complex    *16 a(lda,*), acopy(n,*), cwk(n,*)
        end subroutine
      end interface

      interface
        subroutine de2bsb (n, mxeval, a, lda, ncoda, elow,ehigh, neval, eval, acopy, wk, iwk)
          integer    n, mxeval, lda, ncoda, neval, iwk(*)
          double precision elow, ehigh, a(lda,*), eval(*),                  &
     &           acopy(ncoda+1,*), wk(n,*)
        end subroutine
      end interface

      interface
        subroutine de2bsf (n, mxeval, a, lda, elow, ehigh,neval, eval, acopy, wk, iwk)
          integer    n, mxeval, lda, neval, iwk(*)
          double precision elow, ehigh, a(lda,*), eval(mxeval),             &
     &           acopy(n,n), wk(n,5)
        end subroutine
      end interface

      interface
        subroutine de2ccg (n, a, lda, eval, evec, ldevec,acopy, rwk, cwk)
          integer    n, lda, ldevec
          double precision rwk(*)
          complex    *16 a(lda,*), eval(*), evec(ldevec,*), acopy(n,*),     &
     &           cwk(n,*)
        end subroutine
      end interface

      interface
        subroutine de2cch (n, a, lda, eval, evec, ldevec,acopy, cwork)
          integer    n, lda, ldevec
          complex    *16 a(lda,*), eval(*), evec(ldevec,*), acopy(n,*),     &
     &           cwork(n,*)
        end subroutine
      end interface

      interface
        subroutine de2chf (n, a, lda, eval, evec, ldevec,acopy, rwk, cwk)
          integer    n, lda, ldevec
          double precision eval(*), rwk(*)
          complex    *16 a(lda,*), evec(ldevec,*), acopy(n,*), cwk(*)
        end subroutine
      end interface

      interface
        subroutine de2crg (n, a, lda, eval, evec, ldevec,acopy, ecopy, rwk)
          integer    n, lda, ldevec
          double precision a(lda,*), acopy(n,*), ecopy(n,*), rwk(n,*)
          complex    *16 eval(*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine de2crh (n, a, lda, eval, evec, ldevec,acopy, ecopy)
          integer    n, lda, ldevec
          double precision a(lda,*), acopy(n,*), ecopy(n,*)
          complex    *16 eval(*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine de2csb (n, a, lda, ncoda, eval, evec,ldevec, acopy, wk)
          integer    n, lda, ncoda, ldevec
          double precision a(lda,*), eval(*), evec(ldevec,*),               &
     &           acopy(ncoda+1,*), wk(*)
        end subroutine
      end interface

      interface
        subroutine de2ehf (n, nevec, a, lda, small, eval,evec, ldevec, acopy, rw1, rw2, cwk, iwk)
          integer    n, nevec, lda, ldevec, iwk(*)
          double precision eval(*), rw1(n,*), rw2(n,*)
          complex    *16 a(lda,*), evec(ldevec,*), acopy(n,*), cwk(n,*)
          logical    small
        end subroutine
      end interface

      interface
        subroutine de2esb (n, nevec, a, lda, ncoda, small,eval, evec, ldevec, acopy, wk, iwk)
          integer    n, nevec, lda, ncoda, ldevec, iwk(*)
          double precision a(lda,*), eval(*), evec(ldevec,*),               &
     &           acopy(ncoda+1,*), wk(*)
          logical    small
        end subroutine
      end interface

      interface
        subroutine de2fhf (n, mxeval, a, lda, elow, ehigh,neval, eval, evec, ldevec, acopy, ecopy, rwk, cwk, iwk)
          integer    n, mxeval, lda, neval, ldevec, iwk(*)
          double precision elow, ehigh, eval(*), ecopy(n,*), rwk(n,*)
          complex    *16 a(lda,*), evec(ldevec,*), acopy(n,*), cwk(n,*)
        end subroutine
      end interface

      interface
        subroutine de2fsb (n, mxeval, a, lda, ncoda, elow,ehigh, neval, eval, evec, ldevec, acopy, wk1, wk2, iwk)
          integer    n, mxeval, lda, ncoda, neval, ldevec, iwk(*)
          double precision elow, ehigh, a(lda,*), eval(*), evec(ldevec,*),  &
     &           acopy(ncoda+1,*), wk1(n,*), wk2(*)
        end subroutine
      end interface

      interface
        subroutine de2fsf (n, mxeval, a, lda, elow, ehigh,neval, eval, evec, ldevec, acopy, wk, iwk)
          integer    n, mxeval, lda, neval, ldevec, iwk(*)
          double precision elow, ehigh, a(lda,*), eval(*), evec(ldevec,*),  &
     &           acopy(n,*), wk(n,*)
        end subroutine
      end interface

      interface
        subroutine de2lcg (n, a, lda, eval, acopy, rwk, cwk)
          integer    n, lda
          double precision rwk(*)
          complex    *16 a(lda,*), eval(*), acopy(n,*), cwk(n,*)
        end subroutine
      end interface

      interface
        subroutine de2lch (n, a, lda, eval, acopy)
          integer    n, lda
          complex    *16 a(lda,*), eval(*), acopy(n,*)
        end subroutine
      end interface

      interface
        subroutine de2lhf (n, a, lda, eval, acopy, rwk, cwk)
          integer    n, lda
          double precision eval(*), rwk(*)
          complex    *16 a(lda,*), acopy(n,*), cwk(*)
        end subroutine
      end interface

      interface
        subroutine de2lrg (n, a, lda, eval, acopy, rwk)
          integer    n, lda
          double precision a(lda,*), acopy(n,*), rwk(n,*)
          complex    *16 eval(*)
        end subroutine
      end interface

      interface
        subroutine de2lrh (n, a, lda, eval, acopy)
          integer    n, lda
          double precision a(lda,*), acopy(n,*)
          complex    *16 eval(*)
        end subroutine
      end interface

      interface
        subroutine de2lsb (n, a, lda, ncoda, eval, acopy, wk)
          integer    n, lda, ncoda
          double precision a(lda,*), eval(*), acopy(ncoda+1,*), wk(*)
        end subroutine
      end interface

      interface
        subroutine de3crg (n, low, igh, a, lda, scale)
          integer    n, low, igh, lda
          double precision a(lda,*), scale(*)
        end subroutine
      end interface

      interface
        subroutine de3crh (n, low, igh, a, lda, eval, ecopy,ldecop, evec, ldevec, vector)
          integer    n, low, igh, lda, ldecop, ldevec
          double precision a(lda,*), ecopy(ldecop,*)
          complex    *16 eval(*), evec(ldevec,*)
          logical    vector
        end subroutine
      end interface

      interface
        subroutine de3lsf (n, a, lda, d, e, e2)
          integer    n, lda
          double precision a(lda,*), d(*), e(*), e2(*)
        end subroutine
      end interface

      interface
        subroutine de4crg (n, low, igh, a, lda, ort, work)
          integer    n, low, igh, lda
          double precision a(lda,*), ort(*), work(*)
        end subroutine
      end interface

      interface
        subroutine de4esf (n, nevec, a, evec, ldevec, e)
          integer    n, nevec, ldevec
          double precision a(n,*), evec(ldevec,*), e(*)
        end subroutine
      end interface

      interface
        subroutine de5crg (n, low, igh, a, lda, ort, evec,ldevec, g)
          integer    n, low, igh, lda, ldevec
          double precision a(lda,*), ort(*), evec(ldevec,*), g(*)
        end subroutine
      end interface

      interface
        subroutine de7crg (n, nevec, eval, evec, ldevec,vector)
          integer    n, nevec, ldevec
          complex    *16 eval(*), evec(ldevec,*)
          logical    vector
        end subroutine
      end interface

      interface
        subroutine dg2ccg (n, a, lda, b, ldb, alpha, beta,evec, ldevec, acopy, bcopy)
          integer    n, lda, ldb, ldevec
          complex    *16 a(lda,*), b(ldb,*), alpha(*), beta(*),             &
     &           evec(ldevec,*), acopy(n,*), bcopy(n,*)
        end subroutine
      end interface

      interface
        subroutine dg2crg (n, a, lda, b, ldb, alpha, beta,evec, ldevec, acopy, bcopy, ecopy, rwk, cwk)
          integer    n, lda, ldb, ldevec
          double precision a(lda,*), b(ldb,*), beta(*), acopy(n,*),         &
     &           bcopy(n,*), ecopy(n,*), rwk(*)
          complex    *16 alpha(*), evec(ldevec,*), cwk(*)
        end subroutine
      end interface

      interface
        subroutine dg2lcg (n, a, lda, b, ldb, alpha, beta,acopy, bcopy)
          integer    n, lda, ldb
          complex    *16 a(lda,*), b(ldb,*), alpha(*), beta(*),             &
     &           acopy(n,*), bcopy(n,*)
        end subroutine
      end interface

      interface
        subroutine dg2lrg (n, a, lda, b, ldb, alpha, beta,acopy, bcopy, rwk, cwk)
          integer    n, lda, ldb
          double precision a(lda,*), b(ldb,*), beta(*), acopy(n,*),         &
     &           bcopy(n,*), rwk(*)
          complex    *16 alpha(*), cwk(*)
        end subroutine
      end interface

      interface
        subroutine dg3ccg (n, a, b, ijob, evec, ldevec)
          integer    n, ijob, ldevec
          complex    *16 a(n,*), b(n,*), evec(ldevec,*)
        end subroutine
      end interface

      interface
        subroutine dg4ccg (n, a, b, evec, ldevec, ijob,alpha, beta)
          integer    n, ldevec, ijob
          complex    *16 a(n,*), b(n,*), evec(ldevec,*), alpha(*), beta(*)
        end subroutine
      end interface

      interface
        subroutine dg5ccg (n, nevec, alpha, beta, evec,ldevec, vector, eval)
          integer    n, nevec, ldevec
          complex    *16 alpha(*), beta(*), evec(ldevec,*), eval(*)
          logical    vector
        end subroutine
      end interface

      interface
        subroutine dg7crg (n, nevec, alpha, beta, evec,ldevec, vector, eval)
          integer    n, nevec, ldevec
          double precision beta(*)
          complex    *16 alpha(*), evec(ldevec,*), eval(*)
          logical    vector
        end subroutine
      end interface

      interface
        subroutine divpbs (ido, neq, fcn, x, xend, tol,param, y)
          integer    ido, neq
          double precision x, xend, tol, param(*), y(*)
          external   fcn
        end subroutine
      end interface

      interface
        subroutine di2pbs (ido, n, fcn, x, xend, tol, param,y, r, s, vnorm, wk)
          integer    ido, n
          double precision x, xend, tol, param(*), y(*), r(*), s(*), wk(*)
          external   fcn, vnorm
        end subroutine
      end interface

end module
