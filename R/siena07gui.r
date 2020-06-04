#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: siena07Gui.r
# *
# * Description: This module contains the code controlling the gui for siena07.
# *
# *****************************************************************************/
##@siena07Gui siena07 Create program screen for siena07
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
        tt <- tktoplevel()
    }
    tkwm.title(tt,'Siena07')
    frame <- tkframe(tt, width=300, height=300, relief='ridge')
    tkpack(frame, side='top', padx=5)
    button1 <- tkbutton(frame, command=myInterrupt, text='Interrupt')
    button2 <- tkbutton(frame, command=myEndPhase2, text='End Phase2',
                       state='disabled')
    button3 <- tkbutton(frame, command=myRestart, text='Restart')
    tkgrid.configure(button1, column=1, columnspan=2, row=1, padx=20, pady=20)
    tkgrid.configure(button2, column=3, row=1, padx=20)
    tkgrid.configure(button3, column=4, row=1, padx=20)
    phaselabel <- tklabel(frame, text='Phase')
    subphaselabel <- tklabel(frame, text='Subphase', state='disabled')
    iterationlabel <- tklabel(frame, text='Iteration')
    label1 <- tklabel(frame, text='ProgressBar')

    phase <- tkentry(frame, width=2)

    subphase <- tkentry(frame, width=2, state='disabled')
    iteration <- tkentry(frame, width=6)
    progressbar <- ttkprogressbar(frame, max=2000, length=120)

    tkgrid.configure(phaselabel, column=1, row=2, pady=5)
    tkgrid.configure(subphaselabel, column=2, row=2)
    tkgrid.configure(iterationlabel, column=3, row=2)
    tkgrid.configure(label1, column=4, row=2, padx=5)
    tkgrid.configure(phase, column=1, row=3, pady=3)
    tkgrid.configure(subphase, column=2, row=3, padx=10)
    tkgrid.configure(iteration, column=3, row=3, padx=10)

    tkgrid.configure(progressbar, column=4, padx=5, row=3)
    label2 <- tklabel(frame, text='Current parameter values')
    label3 <- tklabel(frame, text='Quasi-autocorrelations')
    label4 <- tklabel(frame, text='Deviation values')

    tkgrid.configure(label2, column=1, columnspan=2, row=4, padx=10)
    tkgrid.configure(label3, column=3, row=4, padx=10)
    tkgrid.configure(label4, column=4, row=4, padx=10)

    text1 <- tktext(frame, height=6, width=14)

    text2 <- tktext(frame, height=6, width=14)
    text3 <- tktext(frame, height=6, width=14)
    tkgrid.configure(text1, column=1, columnspan=2, row=5, padx=20, pady=5)

    tkgrid.configure(text2, column=3, row=5, padx=20)
    tkgrid.configure(text3, column=4, row=5, padx=20)
    ilcampo <- tclVar()
    tcl("image", "create", "photo", ilcampo, file=imagepath)
    frame2 <- tkframe(tt, width=300, height=300, relief='ridge')
    tkpack(frame2, side='bottom', padx=5)
    imgAsLabel <- tklabel(frame2, image=ilcampo)
    tkgrid.configure(imgAsLabel, pady=10)
    tkinsert(phase, 0, ' 1')
    tkgrab.set(tt)
    tcl('update')
    tkfocus(tt)
    list(tt=tt, pb=progressbar, earlyEndPhase2=button2, current=text1,
         quasi=text2, deviations=text3, phase=phase, subphase=subphase,
         iteration=iteration, subphaselabel=subphaselabel,
		 phaselabel=phaselabel)
}


#tkconfigure(button2,state='normal')


