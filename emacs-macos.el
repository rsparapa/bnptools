
; copy these settings to your ~/.emacs file

; by default,the Command (Apple or Mac) key is meta and Option is alt
; reverse the setting so Cut/Copy/Paste/etc. are like other Mac software
; that follows the Apple Human Interface Guidelines

; the modifier key assumed when the Command key is pressed
(setq mac-command-modifier 'alt)

; the modifier key assumed when the Option key is pressed
(setq mac-option-modifier 'meta)

(defun open-directory ()
  "*Open a directory in the Finder."
  (interactive)
  (call-process "open" nil nil nil
		(read-string "Open directory:" default-directory))
)

(defun clipboard-kill-region-if-active (&optional BEG END) 
"*Cut: workaround for bug/feature."
  (interactive)
  (if (use-region-p) (clipboard-kill-region (region-beginning) (region-end)))
)

; Mac command keys using alt modifier
; as of emacs 27/28, you need to delete the key definition before setting
(define-key key-translation-map [(alt ?x)] nil) 
(global-set-key [(alt x)] 'clipboard-kill-region-if-active)  ; cut
(define-key key-translation-map [(alt ?c)] nil)
(global-set-key [(alt c)] 'clipboard-kill-ring-save)         ; copy
(define-key key-translation-map [(alt ?v)] nil)
(global-set-key [(alt v)] 'clipboard-yank)                   ; paste
(define-key key-translation-map [(alt ?a)] nil)
(global-set-key [(alt a)] 'mark-whole-buffer)                ; all
(define-key key-translation-map [(alt ?k)] nil)
(global-set-key [(alt k)] 'kill-this-buffer)                 ; kill
(define-key key-translation-map [(alt ?o)] nil)
(global-set-key [(alt o)] 'open-directory)                   ; open
(define-key key-translation-map [(alt ?q)] nil)
(global-set-key [(alt q)] 'save-buffers-kill-terminal)       ; quit
(define-key key-translation-map [(alt ?s)] nil)
(global-set-key [(alt s)] 'save-buffer)                      ; save
