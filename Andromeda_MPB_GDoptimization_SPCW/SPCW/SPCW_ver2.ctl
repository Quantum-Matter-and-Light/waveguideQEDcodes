;; This CTL file is to calculate the dispersion relation for "Si3N4 Squircle W1 waveguide structure" by MPB
;; Structure parameters : (1) Lattice constant : a
;;			  (2) Air-hole radius (r1: Squircle size) (r2: Second-line holes) (r3: all the other holes)
;;			  (3) Structure thickness : thk
;;			  (4) Width of air-line-slot: wid
;;			  (5) Ellipcity of squircles : Alpha (the larger value you choose, the more elongated squircle in x-direction)
;;			  (6) Shift of the squircles in y-direction : s1 
;;			  (7) Shift of the second line of holes in y-direction : s2
;;			  (8) offset : Additional dielectric slap between the first line squircles and air line-slot 
;;			  (9) Squareness : Nsq
;;			  (*) Sinusoidal modulation of the edge of th air line-slot (Alligator structrue) : Amplitude

;; W1 waveguide structure parameters ;;
(define pi (* 2 (acos 0)))	; pi
(define ix_SiN 1.9935)		; index of Si3N4
(define myvacuum 1)
;(define-param kst 0.42)
;(define-param ked 0.47)
(define-param kstart 0.4)
(define-param kend 0.4997)
;; PhC dimensions ;;
(define-param Nsq 4)  ; (Only even number) Nsq = 2 -> Circle // Nsq = 4 -> Squircle // Nsq =6, 8, 10 -> Square
(define-param amp 0)
(define-param wid 114)
(define-param thk 200)
(define-param a 373)
(define-param y 8) ; the number of unitcell in y-direction
(define-param r1 110) ; radius of the 1st air-holes
(define-param r2 110) ; radius of the 2nd air-holes
(define-param r3 110) ; radius of all the other air-holes
(define-param s1 0) ;minimum distance between air hole and air slot
(define-param s2 0) ; minimum distance between 1st air hole and 2nd air hole
(define-param offset 100)
(define-param slot (* 2 (+ amp wid)))
(define-param alpha 1.15) ; ellipcity of the 2nd holes

;; define unit conversion function ;;
(define (simUnits x) (/ x a))

;; Unit conversion ;;
(define-param MPBthk (simUnits thk))
(define-param MPBa (simUnits a))
(define-param MPBr1 (simUnits r1))
(define-param MPBr2 (simUnits r2))
(define-param MPBr3 (simUnits r3))
(define-param MPBs1 (simUnits s1))
(define-param MPBs2 (simUnits s2))
(define-param MPBslot (simUnits slot))
(define-param MPBoffset (simUnits offset))

(define-param elpta (* alpha r1))
(define-param elptb (/ r1 alpha))
(define-param MPBelpta (simUnits elpta))
(define-param MPBelptb (simUnits elptb))

(define-param wX 1) ; width of blocks in x-direction
(define-param L1 (+ MPBs1 MPBoffset (/ (sqrt 3) 2))) ; distance between central line and the center of 1st air hole (squircle)
(define-param L2 (+ MPBs2 MPBoffset (sqrt 3))) ; distance between central line and the center of 2nd air hole
(define-param wY (+ (* (sqrt 3) y) (* 2 L2))) ; width of blocks in y-direction
(define-param sy wY) ; the computational size in y-axis

;; AUX Setting ;;
;; AUX settings ;;
(define-param res 20)
(set-param! resolution res)
(define dx (/ 1 res 10))
(define-param kpts 10)
(define-param nbands 24)

;; Cell dimensions ;;
(define sx wX)
(define-param sz 10)
(set! geometry-lattice (make lattice (size sx sy sz)))

;; Create geometry ;;
(define (drawBlock x)
  (let ((dy (+ wid (* amp (cos (* 2 pi x))))))
       (let ((cy (/ dy 2)))
            (list
              (make block
                (center x (simUnits cy) 0)
                (size (* 1.1 dx) (simUnits dy) (simUnits thk))
                (material (make dielectric (index myvacuum)))
              )
              (make block
                (center x (* -1 (simUnits cy)) 0)
                (size (* 1.1 dx) (simUnits dy) (simUnits thk))
                (material (make dielectric (index myvacuum)))
              )
            )
        )
))

(define (drawSquircle x)
  (let ((dy (* 2 MPBelptb (expt (- 1 (/ (expt x 4) (expt MPBelpta 4))) (/ 1 4)))))
       (list 
	 (make block
	  (center (- x 0.5) L1 0)
	  (size (* 1.1 dx) dy MPBthk)
	  (material (make dielectric (index myvacuum)))
	 )
	 (make block
	  (center (- x 0.5) (* -1 L1) 0)
	  (size (* 1.1 dx) dy MPBthk)
	  (material (make dielectric (index myvacuum)))
	 )
	 (make block
	  (center (+ x 0.5) L1 0)
	  (size (* 1.1 dx) dy MPBthk)
	  (material (make dielectric (index myvacuum)))
	 )
	 (make block
	  (center (+ x 0.5) (* -1 L1) 0)
	  (size (* 1.1 dx) dy MPBthk)
	  (material (make dielectric (index myvacuum)))
	 )
       )
)) 

(set! geometry '())
(set! geometry
 (append
  (list
   (make block (center 0 0 0) (size sx sy MPBthk) (material (make dielectric (index ix_SiN)))) ; dielectric Block
   ;(make cylinder (center 0.5 L1 0) (radius MPBr1) (height MPBthk) (material (make dielectric (index myvacuum))))  ; First air hole right above
   ;(make cylinder (center -0.5 L1 0) (radius MPBr1) (height MPBthk) (material (make dielectric (index myvacuum)))) ; First air hole left above
   ;(make cylinder (center 0.5 (* -1 L1) 0) (radius MPBr1) (height MPBthk) (material (make dielectric (index myvacuum)))) ; First air hole right below
   ;(make cylinder (center -0.5 (* -1 L1) 0) (radius MPBr1) (height MPBthk) (material (make dielectric (index myvacuum)))) ; First air hole left below
   (make cylinder (center 0 L2 0) (radius MPBr2) (height MPBthk) (material (make dielectric (index myvacuum))))  ; First air hole right above
   (make cylinder (center 0 (* -1 L2) 0) (radius MPBr2) (height MPBthk) (material (make dielectric (epsilon 1)))) ; First air hole left above
  )
  (geometric-object-duplicates (vector3 0 (sqrt 3) 0) 0 (/ y 2)
   (make cylinder (center 0.5 (+ (/ (sqrt 3) 2) L2) 0) (radius MPBr3) (height MPBthk)
    (material (make dielectric (index myvacuum))))
  )
  (geometric-object-duplicates (vector3 0 (sqrt 3) 0) 0 (/ y 2)
   (make cylinder (center -0.5 (+ (/ (sqrt 3) 2) L2) 0) (radius MPBr2) (height MPBthk)
    (material (make dielectric (index myvacuum))))
  )
  (geometric-object-duplicates (vector3 0 (sqrt 3) 0) 0 (/ y 2)
   (make cylinder (center 0 (+ (sqrt 3) L2) 0) (radius MPBr2) (height MPBthk)
    (material (make dielectric (index myvacuum))))
  )
  (geometric-object-duplicates (vector3 0 (* -1 (sqrt 3)) 0) 0 (/ y 2)
   (make cylinder (center 0.5 (- (/ (sqrt 3) -2) L2) 0) (radius MPBr2) (height MPBthk)
    (material (make dielectric (index myvacuum))))
  )
  (geometric-object-duplicates (vector3 0 (* -1 (sqrt 3)) 0) 0 (/ y 2)
   (make cylinder (center -0.5 (- (/ (sqrt 3) -2) L2) 0) (radius MPBr2) (height MPBthk)
    (material (make dielectric (index myvacuum))))
  )
  (geometric-object-duplicates (vector3 0 (* -1 (sqrt 3)) 0) 0 (/ y 2)
   (make cylinder (center 0 (- (* -1 (sqrt 3)) L2) 0) (radius MPBr2) (height MPBthk)
    (material (make dielectric (index myvacuum)))) ; First middle hole below
  )  
))

(define (makeModule x)
  (if (>= x 1)
        '()
      (begin
        (set! geometry (append geometry (drawBlock x)))
	(makeModule (+ x dx))    			
      )
  )
)
(define (makeSquircle x)
  (if (>= x MPBelpta)
        '()
      (begin
	(set! geometry (append geometry (drawSquircle x)))
	(makeSquircle (+ x dx))
      )
  )
)				 

(makeModule 0)
(makeSquircle (* -1 MPBelpta))

(define pork (vector3 0.4998 0 0))
(define chicken (vector3 0.4999 0 0))
(define beef (vector3 0.5 0 0))
(define Delicacy (list pork chicken beef))
(define Gamma (vector3 kstart 0 0))
(define M (vector3 kend 0 0))
(define Normal (interpolate kpts (list Gamma M)))
(set! k-points (append Normal Delicacy))
(set! num-bands nbands)
(set! tolerance 1e-9)
;(run-tm display-group-velocities fix-efield-phase output-efield)
(run-zeven display-group-velocities (output-at-kpoint beef fix-efield-phase output-efield-x output-efield-y output-efield-z ))

