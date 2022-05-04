;; Utilities I want to share between exercises.

;; Printing

(load "exdisplay.scm")

(define (replace-all haystack needle replacement)
  (let ((haystack (string->list haystack))
        (replacement (reverse
                      (string->list replacement)))
        (needle-len (string-length needle)))
    (let loop ((haystack haystack) (acc '()))
      (cond ((null? haystack)
             (list->string (reverse acc)))

            ((string-prefix? needle (list->string haystack))
             (loop (list-tail haystack needle-len)
                   (append replacement acc)))

            (else
             (loop (cdr haystack) (cons (car haystack) acc)))))))

;; Generates a properly formatted string of LaTeX.
(define (->tex* expr)
  (let* ((tex-string (expression->tex-string
                      ((prepare-for-printing expr simplify))))
         (len (string-length tex-string)))
    (substring tex-string 10 (- len 3))))

;; Prints the string as a LaTeX code block.
(define (->write-tex tex-string)
  (write-string
   (string-append "\\[ " tex-string " \\]")))

;; Prints the TeX representation of the supplied expression to the screen.
(define (->tex expr)
  (->write-tex (->tex* expr)))

;; Prints an equation code block containing the expression as LaTeX.
(define (->tex-equation* expr #!optional label)
  (string-append
   "\\begin{equation}\n"
   (->tex* expr)
   (if (default-object? label)
       ""
       (string-append "\n\\label{" label "}"))
   "\n\\end{equation}"))

(define (->tex-equation expr #!optional label)
  (write-string (->tex-equation* expr label)))
