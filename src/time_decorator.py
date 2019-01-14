import time

# def timer(start_time,str="Function"):
def timer(func,start_time,str,*args,**kwargs):
    """This function is running the main function ``func`` and times its process.
    
    :param func: any function
    :param start_time: starting time of the main program as given by main program
    :param \*args: arguments for func function
    :param \*\*kwargs: keyword arguments for func function

    :rtype: return of func if func returns anything, else nothing

    
    Note to myself: that func.__name__ returns the name of func as a string
    This could be useful later

    """

    # Start time of sub process
    time_i = time.time()

    # Print Start function statement
    print()
    print("Start %s..." % str)

    # ----- Run function ------
    if callable(func):
        result = func(*args,**kwargs)
    else:
        func(*args,**kwargs)
    # -------------------------

    # End time of sub process
    time_f = time.time()

    # time difference in subprocess
    # Total time needed:
    dtime = time_f - time_i
    
    # Print Time from start of program and Time eeded for Meshing 
    print("Finished %s." % str)
    print("Time: %s sec -- dt: %s sec" % ((time.time() - start_time) , dtime))
    print()

    if callable(func):
        return result
    else:
        return