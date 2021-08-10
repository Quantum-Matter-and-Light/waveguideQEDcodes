;; 2017.08.23 Youn Seok Lee
;; This CTL file is to calculate Casimir-Polder potential for the ground state of Cesium atom near the "Si3N4 Squircle W1 waveguide structure" by MEEP. "filename_xyr" is for the calculation of CP potential in xy-plane at z=0, and "filename_zr" is for the z-direction at x=y=0. 

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

;; distance units in (um) and unity speed of light (c=1) ;;

;(include "/usr/share/meep/casimir.scm")
(include "/home/uqml/apps/meep-1.3/share/meep/parallel.scm")

;; constant;;
(define pi (* 2 (acos 0)))	; pi
(define Si3N4 1.9935)		; index of Si3N4
(define-param myvacuum 1)
(define-param Sigma 1)

;; PhC dimensions (um);;
(define-param Nsq 4)  	; (Only even number) Nsq = 2 -> Circle // Nsq = 4 -> Squircle // Nsq =6, 8, 10 -> More square-like shape
(define-param amp 0)
(define-param wid 114)
(define-param thk 200)
(define-param a 373)	 ; lattice constant
(define-param y 8) 	 ; the number of unitcell in y-direction
(define-param r1 110)	 ; radius of the 1st air-holes (Squircles)
(define-param r2 110)	 ; radius of the 2nd air-holes 
(define-param r3 110) 	 ; radius of all the other air-holes
(define-param s1 -20) 	 ; minimum distance between air hole and air slot
(define-param s2 -20) 	 ; minimum distance between 1st air hole and 2nd air hole
(define-param offset 120)
(define-param slot (* 2 (+ amp wid)))
(define-param alpha 1.15) 	 ; ellipcity of the 2nd holes
(define-param buffer 500) ; distance between the edge of structure and PML layer

;; define unit conversion function ;;
(define (simUnits x) (/ x a))   ; units of length in (um) 

;; Unit conversion ;;
(define-param MEEPthk (simUnits thk))
(define-param MEEPa (simUnits a))
(define-param MEEPr1 (simUnits r1))
(define-param MEEPr2 (simUnits r2))
(define-param MEEPr3 (simUnits r3))
(define-param MEEPs1 (simUnits s1))
(define-param MEEPs2 (simUnits s2))
(define-param MEEPslot (simUnits slot))
(define-param MEEPoffset (simUnits offset))
(define-param MEEPbuffer (simUnits buffer))
(define-param elpta (* alpha r1))
(define-param elptb (/ r1 alpha))
(define-param MEEPelpta (simUnits elpta))
(define-param MEEPelptb (simUnits elptb))

;; Define perfecty-matched-layer thickness ;;
(define-param dpml 1)

;; Define computational cell dimension ;;
(define supercell-x MEEPa) ; width of blocks in x-direction
(define L1 (+ MEEPs1 MEEPoffset (* (/ (sqrt 3) 2) MEEPa))) ; distance between central line and the center of 1st air hole (squircle)
(define L2 (+ MEEPs2 MEEPoffset (* (sqrt 3) MEEPa))) ; distance between central line and the center of 2nd air hole
(define wY (+ (* (sqrt 3) MEEPa y) (* 2 L2))) ; width of blocks in y-direction
(define supercell-y (+ wY (* 2 MEEPbuffer) (* 2 dpml))) ; the computational size in y-axis
(define supercell-z (+ MEEPthk (* 2 MEEPbuffer) (* 2 dpml)))

;; AUX Settings ;;
(define-param res 32)
(define-param crnt 0.5)
(set-param! resolution res)
(set-param! Courant crnt)
(define-param endtime 10)
(define-param kpt (* 0.5 (/ 1 MEEPa)))
(define dx (/ 1 res 10))

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
  (let ((dy (* 2 MEEPelptb (expt (- 1 (/ (expt x Nsq) (expt MEEPelpta 4))) (/ 1 Nsq)))))
       (list 
	 (make block
	  (center x L1 0)
	  (size (* 1.1 dx) dy MEEPthk)
	  (material (make dielectric (index myvacuum)))
	 )
	 (make block
	  (center x (* -1 L1) 0)
	  (size (* 1.1 dx) dy MEEPthk)
	  (material (make dielectric (index myvacuum)))
	 )
       )
)) 

(define (makeModule x)
  			(if (>= x 1)
        		'()
      			(begin (set! geometry (append geometry (drawBlock x)))
			       (makeModule (+ x dx))    			
      			)
  		))
(define (makeSquircle x)
  			(if (>= x MEEPelpta)
       			'()
      			(begin (set! geometry (append geometry (drawSquircle x)))
			       (makeSquircle (+ x dx))
      			)
  		))


;; Sources ;;
(define-param fcen (/ 1 0.852))			; unit in (c/a=1/(lambda/a))
(define-param df 1)

;; position of interest ;;
(define-param xbgn 0)
(define-param ybgn 0)
(define-param zbgn 0)
(define-param r-min -0.2)
(define-param r-max (+ -0.2 (/ 20 res)))
(define-param dr (/ 1 res))
(define r-list (parallel-make-list r-min r-max dr))
(define eps-list (list Si3N4 myvacuum))

(define pol-list (list Ex Ey Ez))

(define Ext-list (list r-list eps-list))
(define Int-list (list pol-list))

(define param-info (make-param-list Ext-list Int-list))
(print param-info)
(define param-list (first param-info))
(define Next (second param-info))
(define Nint (third param-info))
(define Nsims (* Next Nint))

(define gamma-list (make-list Next 0))
(print "Number of simulations: "Nsims"\n")

;%%%%%%%%%%%%%%%%%%% LOCAL COMMUNICATION

;(define nproc (meep-count-processors))
;(define ngroups (min Nsims nproc))
;(define mygroup (meep-divide-parallel-processes ngroups))

;(print "nproc = "nproc"\n")
;(print "ngroups = "ngroups"\n")
;(print "mygroup = "mygroup"\n")

;(define my-sims (get-indices Nsims nproc mygroup)) ;a list of simulations for the group
;(print "Total: my-sims = "my-sims"\n")

(define polstring "xyddz")
(define xstring "")
(define rstring "")
(define prev-r r-min)

(define (run-sim current-sim)
	(let* (  (index-info (get-ie-indices current-sim Next Nint)) 
	 	 (i-internal (first index-info))
	 	 (i-external (second index-info))
		 ;%%%%% Get current simulation parameters
		 (curr-params (list-ref param-list current-sim))
		 (curr-r (first curr-params))
		 (curr-eps (second curr-params))
		 (curr-pol (third curr-params))
		 (ft (meep-type curr-pol))
		 (dt (/ Courant resolution))
	      )
	(print "index-info "index-info", i-internal "i-internal", i-external "i-external", curr-params "curr-params"\n")

		(set! pml-layers (list
                		(make pml (direction Z) (thickness dpml))
		                (make pml (direction Y) (thickness dpml))
		                ;(make pml (direction X) (thickness dpml))
		         )
		)
		
		(set! ensure-periodicity true)
		(set! geometry-lattice (make lattice (size supercell-x supercell-y supercell-z)))
		;Bloch perioodicity
		(set-param! k-point (vector3 kpt 0 0))

		(set! geometry '())

		(set! geometry
		 (append
		  (list
		   (make block (center 0 0 0) (size infinity infinity infinity) (material (make dielectric (index myvacuum))))
 		   (make block (center 0 0 0) (size MEEPa wY MEEPthk) (material (make dielectric (index curr-eps)))) ; dielectric Block   		  
  		   (make cylinder (center (* 0.5 MEEPa) L2 0) (radius MEEPr2) (height MEEPthk) (material (make dielectric (index myvacuum))))  ; First air hole right above
  		   (make cylinder (center (* -0.5 MEEPa) L2 0) (radius MEEPr2) (height MEEPthk) (material (make dielectric (index myvacuum))))  ; First air hole right above
 		   (make cylinder (center (* 0.5 MEEPa) (* -1 L2) 0) (radius MEEPr2) (height MEEPthk) (material (make dielectric (index myvacuum)))) ; First air hole left above
  		   (make cylinder (center (* -0.5 MEEPa) (* -1 L2) 0) (radius MEEPr2) (height MEEPthk) (material (make dielectric (index myvacuum)))) ; First air hole left above
  		  )
  		  (geometric-object-duplicates (vector3 0 (* (sqrt 3) MEEPa) 0) 0 (/ (- y 2) 2)
   			(make cylinder (center 0 (+ (* (/ (sqrt 3) 2) MEEPa) L2) 0) (radius MEEPr3) (height MEEPthk)
   			(material (make dielectric (index myvacuum))))
 		  )
  		  (geometric-object-duplicates (vector3 0 (* -1 (* (sqrt 3) MEEPa)) 0) 0 (/ (- y 2) 2)
  			(make cylinder (center 0 (- (* (/ (sqrt 3) -2) MEEPa) L2) 0) (radius MEEPr2) (height MEEPthk)
    			(material (make dielectric (index myvacuum))))
  		  )
  		  (geometric-object-duplicates (vector3 0 (* (sqrt 3) MEEPa) 0) 0 (/ (- y 2) 2)
   			(make cylinder (center (* 0.5 MEEPa) (+ (* (sqrt 3) MEEPa) L2) 0) (radius MEEPr2) (height MEEPthk)
    			(material (make dielectric (index myvacuum))))
  		  )
  		  (geometric-object-duplicates (vector3 0 (* (sqrt 3) MEEPa) 0) 0 (/ (- y 2) 2)
   			(make cylinder (center (* -0.5 MEEPa) (+ (* (sqrt 3) MEEPa) L2) 0) (radius MEEPr2) (height MEEPthk)
    			(material (make dielectric (index myvacuum))))
  		  )
  		  (geometric-object-duplicates (vector3 0 (* -1 (* (sqrt 3) MEEPa)) 0) 0 (/ (- y 2) 2)
   			(make cylinder (center (* 0.5 MEEPa) (- (* -1 (* (sqrt 3) MEEPa)) L2) 0) (radius MEEPr2) (height MEEPthk)
    			(material (make dielectric (index myvacuum)))) ; First middle hole below
  		  )  
  		  (geometric-object-duplicates (vector3 0 (* -1 (* (sqrt 3) MEEPa)) 0) 0 (/ (- y 2) 2)
   			(make cylinder (center (* -0.5 MEEPa) (- (* -1 (* (sqrt 3) MEEPa)) L2) 0) (radius MEEPr2) (height MEEPthk)
    			(material (make dielectric (index myvacuum)))) ; First middle hole below
  		  )
		))
		
		
		(makeModule 0)
		(makeSquircle (* -1 MEEPelpta))

		
		(if (= ft E-stuff) (begin 
					(set! global-D-conductivity Sigma)
					(set! global-B-conductivity 0)))


		(set! sources (list (make source ( src 
							(make custom-src 
							(src-func (lambda (t) (/ 1 dt)))
							(start-time (* -.25 dt))
							(end-time (* .75 dt))
							(width dt)
							(is-integrated? false))
							;(make gaussian-src (frequency fcen) (fwidth df))
							
						 )
						 (component curr-pol)
						 (center xbgn curr-r 0)
						 (size 0 0 0))
					)
		)

(print "Current polarization: "curr-pol"\n")
(print (string (string-ref polstring curr-pol))"\n")
(print "Current dielectric constant: "curr-eps"\n")
(print "Current position: "curr-r"\n")

	(reset-meep)
	(init-fields)
	(if (= curr-eps myvacuum) (use-output-directory (string-append "out_vac_xyr_xp" (number->string xbgn) "_kx" (number->string kpt))) (use-output-directory (string-append "out_xyr_xp" (number->string xbgn) "_kx" (number->string kpt))))
	(if (not (= curr-r prev-r)) (begin (set! rstring (string-append rstring "r")) (set! prev-r curr-r)))
;	(if (not (= xbgn 0)) (begin (set! xstring (string-append xstring "x")) (set! xbgn 0)))
;	(get-filename-prefix "kpt")
        (run-until endtime
		(at-beginning output-epsilon)
                (after-time 0 (to-appended (string-append (string #\e #\x) (string (string-ref polstring curr-pol)) (string #\k #\x) (number->string kpt) (string #\r) rstring) (at-every dt (in-point (vector3 xbgn curr-r 0) output-efield-x))))
                (after-time 0 (to-appended (string-append (string #\e #\y) (string (string-ref polstring curr-pol)) (string #\k #\x) (number->string kpt) (string #\r) rstring) (at-every dt (in-point (vector3 xbgn curr-r 0) output-efield-y))))
                (after-time 0 (to-appended (string-append (string #\e #\z) (string (string-ref polstring curr-pol)) (string #\k #\x) (number->string kpt) (string #\r) rstring) (at-every dt (in-point (vector3 xbgn curr-r 0) output-efield-z))))
		  )
	)
)

(do ((j 0 (1+ j))) ((= j Nsims)) (run-sim j))