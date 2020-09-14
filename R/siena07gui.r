#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: siena07Gui.r
# *
# * Description: This module contains the code controlling the gui for siena07.
# *
# *****************************************************************************/
##@siena07Gui siena07 Create program screen for siena07
## This is called only if !batch.
siena07Gui <- function(tt, getDocumentation=FALSE)
{
    ##@myInterrupt internal siena07Gui
    myInterrupt<- function()
    {
        UserInterrupt(TRUE)
    }
    ##@myEndPhase2 internal siena07Gui
    myEndPhase2<- function()
    {
        EarlyEndPhase2(TRUE)
    }
    ##@myRestart internal siena07Gui
    myRestart<- function()
    {
        UserRestart(TRUE)
    }
    if (getDocumentation)
    {
        return(getInternals())
    }
    if (is.null(tt))
    {
        ## require(tcltk)
        tt <- tcltk::tktoplevel()
    }
    tcltk::tkwm.title(tt,'Siena07')
    frame <- tcltk::tkframe(tt, width=300, height=300, relief='ridge')
    tcltk::tkpack(frame, side='top', padx=5)
    button1 <- tcltk::tkbutton(frame, command=myInterrupt, text='Interrupt')
    button2 <- tcltk::tkbutton(frame, command=myEndPhase2, text='End Phase2',
                       state='disabled')
    button3 <- tcltk::tkbutton(frame, command=myRestart, text='Restart')
    tcltk::tkgrid.configure(button1, column=1, columnspan=2, row=1, padx=20, pady=20)
    tcltk::tkgrid.configure(button2, column=3, row=1, padx=20)
    tcltk::tkgrid.configure(button3, column=4, row=1, padx=20)
    phaselabel <- tcltk::tklabel(frame, text='Phase')
    subphaselabel <- tcltk::tklabel(frame, text='Subphase', state='disabled')
    iterationlabel <- tcltk::tklabel(frame, text='Iteration')
    label1 <- tcltk::tklabel(frame, text='ProgressBar')

    phase <- tcltk::tkentry(frame, width=2)

    subphase <- tcltk::tkentry(frame, width=2, state='disabled')
    iteration <- tcltk::tkentry(frame, width=6)
    progressbar <- tcltk::ttkprogressbar(frame, max=2000, length=120)

    tcltk::tkgrid.configure(phaselabel, column=1, row=2, pady=5)
    tcltk::tkgrid.configure(subphaselabel, column=2, row=2)
    tcltk::tkgrid.configure(iterationlabel, column=3, row=2)
    tcltk::tkgrid.configure(label1, column=4, row=2, padx=5)
    tcltk::tkgrid.configure(phase, column=1, row=3, pady=3)
    tcltk::tkgrid.configure(subphase, column=2, row=3, padx=10)
    tcltk::tkgrid.configure(iteration, column=3, row=3, padx=10)

    tcltk::tkgrid.configure(progressbar, column=4, padx=5, row=3)
    label2 <- tcltk::tklabel(frame, text='Current parameter values')
    label3 <- tcltk::tklabel(frame, text='Quasi-autocorrelations')
    label4 <- tcltk::tklabel(frame, text='Deviation values')

    tcltk::tkgrid.configure(label2, column=1, columnspan=2, row=4, padx=10)
    tcltk::tkgrid.configure(label3, column=3, row=4, padx=10)
    tcltk::tkgrid.configure(label4, column=4, row=4, padx=10)

    text1 <- tcltk::tktext(frame, height=6, width=14)

    text2 <- tcltk::tktext(frame, height=6, width=14)
    text3 <- tcltk::tktext(frame, height=6, width=14)
    tcltk::tkgrid.configure(text1, column=1, columnspan=2, row=5, padx=20, pady=5)

    tcltk::tkgrid.configure(text2, column=3, row=5, padx=20)
    tcltk::tkgrid.configure(text3, column=4, row=5, padx=20)
    ilcampo <- tcltk::tclVar()
    tcltk::tcl("image", "create", "photo", ilcampo, file=imagepath)
    frame2 <- tcltk::tkframe(tt, width=300, height=300, relief='ridge')
    tcltk::tkpack(frame2, side='bottom', padx=5)
    imgAsLabel <- tcltk::tklabel(frame2, image=ilcampo)
    tcltk::tkgrid.configure(imgAsLabel, pady=10)
    tcltk::tkinsert(phase, 0, ' 1')
    tcltk::tkgrab.set(tt)
    tcltk::tcl('update')
    tcltk::tkfocus(tt)
    list(tt=tt, pb=progressbar, earlyEndPhase2=button2, current=text1,
         quasi=text2, deviations=text3, phase=phase, subphase=subphase,
         iteration=iteration, subphaselabel=subphaselabel,
		 phaselabel=phaselabel)
}


#tcltk::tkconfigure(button2,state='normal')


