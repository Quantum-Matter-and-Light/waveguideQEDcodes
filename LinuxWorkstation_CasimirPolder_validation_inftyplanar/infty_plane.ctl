;;;;;2017.07.15 Casimir-force calculation in nanofiber trap configuration;;;;;
;;;;;distance units in (um) and unity speed of light (c=1);;;;

;(include "/usr/share/meep/casimir.scm")
(include "/usr/share/meep/parallel.scm")

;;constant;;
(define pi (* 2 (acos 0))) ; pi
(define Si3N4 1.45)     ; index of Si3N4
(define-param myvacuum 1)
(define-param Sigma 1)
(define-param crnt 0.5)
;; Define perfecty-matched-layer thickness ;;
(define-param dpml 1)

;;define unit conversion function;;
(define (simUnits x) (/ x a))

(define supercell-x 8)		; computational cell size in x
(define supercell-y 4)	; computational cell size in y
(define supercell-z 4)	; computational cell size in z

;; AUX Settings ;;
(define-param res 50)
(set-param! resolution res)
(set-param! Courant crnt)
(define-param endtime 10)
(define-param kpt 0)

;; Sources ;;
(define-param fcen (/ 1 0.852))			; unit in (c/a=1/(lambda/a))
(define-param df 1)

;; position of interest ;;
(define-param xbgn 0)
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
		 (currindex (second curr-params))
		 (curr-pol (third curr-params))
		 (ft (meep-type curr-pol))
		 (dt (/ Courant resolution))
		 
	      )
	(print "index-info "index-info", i-internal "i-internal", i-external "i-external", curr-params "curr-params"\n")

		(set! pml-layers (list
                		;(make pml (direction Z) (thickness dpml))
		                ;(make pml (direction Y) (thickness dpml))
		                (make pml (direction X) (thickness dpml))
		         )
		)
		
		(set! ensure-periodicity true)
		(set! geometry-lattice (make lattice (size supercell-x supercell-y supercell-z)))
		;Bloch perioodicity
		(set-param! k-point (vector3 kpt 0 0))
		(set! geometry '()) 	; initialize empty geometry
		(set! geometry (append (list 
					(make block (material (make dielectric (index myvacuum))) (center 0 0 0) (size infinity infinity infinity))
					(make block (material (make dielectric (index currindex))) (center -1.5 0 0) (size 3 infinity infinity)))))

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
						 (center curr-r 0 0)
						 (size 0 0 0))
					)
		)

(print "Current polarization: "curr-pol"\n")
(print (string (string-ref polstring curr-pol))"\n")
(print "Current dielectric constant: "currindex"\n")
(print "Current position: "curr-r"\n")

	(reset-meep)
	(init-fields)
	(if (= currindex myvacuum) (use-output-directory (string-append "out_vac_xyr_" (number->string xbgn))) (use-output-directory (string-append "out_xyr_" (number->string xbgn))))
	(if (not (= curr-r prev-r)) (begin (set! rstring (string-append rstring "r")) (set! prev-r curr-r)))
;	(if (not (= xbgn 0)) (begin (set! xstring (string-append xstring "x")) (set! xbgn 0)))
;	(get-filename-prefix "kpt")
        (run-until endtime
		(at-beginning output-epsilon)
                (after-time 0 (to-appended (string-append (string #\e #\x) (string (string-ref polstring curr-pol)) (string #\k #\x) (number->string kpt) (string #\r) rstring) (at-every dt (in-point (vector3 curr-r 0 0) output-efield-x))))
                (after-time 0 (to-appended (string-append (string #\e #\y) (string (string-ref polstring curr-pol)) (string #\k #\x) (number->string kpt) (string #\r) rstring) (at-every dt (in-point (vector3 curr-r 0 0) output-efield-y))))
                (after-time 0 (to-appended (string-append (string #\e #\z) (string (string-ref polstring curr-pol)) (string #\k #\x) (number->string kpt) (string #\r) rstring) (at-every dt (in-point (vector3 curr-r 0 0) output-efield-z))))
		  )
	)
)

(do ((j 0 (1+ j))) ((= j Nsims)) (run-sim j))
