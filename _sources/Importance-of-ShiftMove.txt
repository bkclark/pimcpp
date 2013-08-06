In PIMC++, time slices are shared among processors like so:
Image:shiftMove1.png

Because time slices at processor boundaries are shared (i.e. both
processors need to use information about that time slice) shared time
slices are prohibited from being moved by any processor. This ensures
that all processors have accurate information about shared a time slice.
If a given slice always stayed between two processors, then it would
never be moved in the simulation and the simulation would be non-ergodic
(and in fact will give pathological results). To ensure this doesn't
happen, there is a special move, "ShiftMove", that shifts the slices so
that shared slices are now in the middle of a processor:
Image:shiftMove2.png

If your algorithm doesn't regularly call ShiftMove, this will never
happen. A general rule is to call ShiftMove every time a BisectionBlock
gets called. The ShiftMove will show up in the time analysis of your
simulation and you can check if you are spending too much time calling
it (more then 1%). In that case, you can always cut back the amount it
is called. **NOTE:** ShiftMove is necessary EVEN when running in serial.
There is still a "stuck" slice between processor 0 and itself.
