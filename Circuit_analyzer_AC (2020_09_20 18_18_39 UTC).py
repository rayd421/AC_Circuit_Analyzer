import numpy as np
import cmath as cm

class comp:
    def __init__(self, name, ctype, value, net1, net2):
        self.name = name
        self.ctype = ctype
        self.value = float(value)
        self.net = [net1, net2]
        self.visited = False
        self.loops = []
        self.polarity = []
        self.voltage = []

    def other_pin(self, nin):
        if self.net[0] == nin:
            return self.net[1]
        else:
            return self.net[0]


class net:
    def __init__(self, name):
        self.name = name
        self.comps = []
        self.voltage =[]


class loop:
    def __init__(self, lid):
        self.id = lid
        self.comps = []
        self.nets = []


comps = []
comp_map = {}
nets = {}
loops = []

def clear_visited():
    global comps
    for c in comps:
        c.visited = False


#
# adds in a class object to the comps
#
def add_class(cl):
    global comps, comp_map
    if cl.name in comp_map:
        print("Hey there, the name has already been used ", cl.name)
        return
    comps.append(cl)
    comp_map[cl.name] = cl


#
# reads in a file
#
def readin(fn):
    with open(fn, "r") as f:
        for l in f:
            if l[0] == "#":
                pass
            elif l[0].upper() == 'R':
                flds = l.strip().split()
                cls = comp(flds[0], "R", flds[1], flds[2], flds[3])
                add_class(cls)
            elif l[0].upper() == 'V':
                flds = l.strip().split()
                cls = comp(flds[0], "V", flds[1], flds[2], flds[3])
                add_class(cls)
            elif l[0].upper() == 'F':
                flds = l.strip().split()
                cls = comp(flds[0], "F", flds[1], None, None)
                add_class(cls)
            elif l[0].upper() == 'C':
                flds = l.strip().split()
                cls = comp(flds[0], "C", flds[1], flds[2], flds[3])
                add_class(cls)
            elif l[0].upper() == 'E':
                break
            else:
                print("I'm so lost, you are messing with me")
                print("---->>>{}<<<----".format(l.strip()))


#
# Make a net list of all the components on a net
#
def makenets():
    global comps, nets
    for c in comps:
        for p in c.net:
            if not (p in nets):
                nets[p] = net(p)
            nets[p].comps.append(c)


def makeloops():
    global comps, nets, loops
    loopnum = 0  # the number of the current loop to make
    for c in comps:  # go through all the components
        if len(c.loops) == 0:  # process any component not in a loop
            print("trying {}".format(c.name))
            clear_visited()  # clear the visited flag in all components
            n = nets[c.net[0]]  # the net to leave this block on (net 0 or + pin)
            c.visited = True  # this block has been visited
            fifo = []  # start with an empty fifo
            found = False  # flag indicating the loop is found
            for c1 in n.comps:  # check all the components on the net from the starting component
                if not (c1.name == c.name):  # check to make sure we don't put the start there twice
                    # the fifo is a tuple of ( component name, net name, component list, net list )
                    fifo.append((c1.name, n.name, [c.name, c1.name], [n.name]))  # put comp in fifo
            while (not found) and len(fifo) > 0:  # run until a loop is found, or the fifo goes empty
                fd = fifo[0]  # get the first entry from the fifo
                fifo = fifo[1:]  # takes front off the fifo
                print("wc is {}".format(fd))
                wc = comp_map[fd[0]]  # get the component (fd[0] is the component name)
                wc.visited = True  # mark this as visited for later use
                other_net = wc.other_pin(fd[1])  # get the net on the other component pin
                print("came in on {}, going out on {}".format(fd[1], other_net))
                onet = nets[other_net]  # get the net class for that net
                for c2 in onet.comps:  # check all the components on the net
                    if c2.name == c.name:  # if back to where we started, we found a loop
                        print("have a path there", fd[2], fd[3])
                        found = True  # so the outer loop will stop unting from the fifo
                        nlist = [x for x in fd[3]]
                        nlist.append(other_net)
                        for ix in range(len(fd[2])):  # place the components on the loop
                            c3 = fd[2][ix]  # get this component name from the list
                            n3 = nlist[ix]  # get the net name from the net list
                            comp_map[c3].loops.append(loopnum)  # added it to the loop
                            if n3 == comp_map[c3].net[0]:  # check the net, and set the polarity
                                comp_map[c3].polarity.append(1)  # on pin 0 is positive
                            else:
                                comp_map[c3].polarity.append(-1)  # on pin 1 is negative
                        lp = loop(loopnum)  # make a new loop
                        lp.comps = [x for x in fd[2]]  # add the components to the loop
                        lp.nets = nlist  # add the nets to the loop
                        loops.append(lp)  # add this loop to the list
                        loopnum = len(loops)  # point to the next loop to add
                        break
                    if not c2.visited:  # Not the start, if not visited, add it to fifo
                        olist = [q for q in fd[2]]  # copy the input component list
                        olist.append(c2.name)  # add this one to the list
                        oname = [q for q in fd[3]]  # copy the signal name list
                        oname.append(onet.name)  # add the signal name to the list
                        fifo.append((c2.name, onet.name, olist, oname))  # put on the fifo


def makeRarray():
    global comps, loops, nets
    nl = len(loops)
    ra = np.array([[complex(0.0) for x in range(nl)] for y in range(nl)])
    print(ra)
    for l in range(nl):
        for c in loops[l].comps:
            co = comp_map[c]
            if co.ctype == 'R':
                for iy in range(len(co.loops)):
                    if iy == l:
                        break
                for ix in range(len(co.loops)):
                    ra[l, co.loops[ix]] += co.value * co.polarity[ix] * co.polarity[iy]
            elif co.ctype == 'C':
                for iy in range(len(co.loops)):
                    if iy == l:
                        break
                for ix in range(len(co.loops)):
                    ra[l, co.loops[ix]] += co.value * co.polarity[ix] * co.polarity[iy]
    return ra


def makeVarray():
    global comps, loops, nets
    nl = len(loops)
    va = np.array([[complex(0.0)] for y in range(nl)])
    for l in range(nl):
        for c in loops[l].comps:
            co = comp_map[c]
            if co.ctype == 'V':
                for ix in range(len(co.loops)):
                    if co.loops[ix] == l:
                        va[l, 0] += co.value * co.polarity[ix]
    return va

def get_volt():
    global comps, ires
    for c in comps:
        a = [0, 0, 0]
        if c.ctype == 'V':
            for l in range(len(c.loops)):               #get voltages of voltage source components
                a[c.loops[l]] = c.value * c.polarity[l]  #store voltage according to loop
        elif c.ctype == 'R' or 'C':                            #get voltage of resistor components
            for k in range(len(c.loops)):               #get voltages of voltage source components
                a[c.loops[k]] = c.value * ires[c.loops[k]]*c.polarity[k]      #store voltage according to loop
        c.voltage = a
        ind = 1
        for li in c.voltage:
            print("The complex voltage across ", c.name, " in loop ", ind, " is ", cm.polar(li))
            ind += 1

def other_term(array, nin):             #function to get other pin net
    if array[0] == nin:
        return array[1]
    else:
        return array[0]


def process_new_net(new_n, com_list, net_list, volt):              #recursive function to process nets
    process_n = new_n
    for i in range(len(nets[process_n].comps)):                    # check and run for new components connected to net
        new_comp = ((nets[process_n]).comps[i])
        if new_comp.name not in com_list:                          # if component is not in already seen list, continue
            new_n = other_term(nets[process_n].comps[i].net, process_n)  # get net of new component
            if new_n in net_list:
                com_list.append(((nets[process_n]).comps[i]).name) # add new component to already processed list
            else:
                volt = volt + sum(new_comp.voltage)                # add voltage of new component
                nets[new_n].voltage = volt                         # store voltage of new net
                print(new_n, " voltage in complex polar form is ", cm.polar(volt))                 # print voltage of net
                com_list.append(((nets[process_n]).comps[i]).name)    # add new component to already processed list
                if new_n not in net_list:                             #if new net not in net list, continue
                    net_list.append(new_n)
                    process_new_net(new_n, com_list, net_list, volt)  #process new net


def get_net_voltages():
    global comps
    com_list = []
    net_list = []
    volt = 0
    for n in comps:                        #Get starting component connected to ground
        if 'gnd' in n.net:
            com_list.append(n.name)         #add starting component to component seen list
            net_list.append('gnd')          #add 'gnd' starting net to net processed list
            new_n = n.other_pin('gnd')      #get other net
            volt = volt + sum(n.voltage)    #add voltage
            print(new_n, " voltage in complex polar form is ", cm.polar(volt))
            break
    process_new_net(new_n, com_list, net_list, volt)



def get_comp_z():       #find the impedance of capacitors
    freq = 0
    for n in comps:
        if n.ctype == 'F':
            freq = n.value
    if freq != 0:
        for c in comps:
            if c.ctype == 'C':
                cap = c.value
                z = (complex(0, ((-1)/((2*(cm.pi))*(freq)*(cap)))))
                c.value = z








readin("simplecir-1.txt")
makenets()
makeloops()
get_comp_z()

# for c in comps:
#   print(c.name,c.value,c.net)
# for n in nets.keys():
#    print(n)
#    for c in nets[n].comps:
#        print("  ",c.name)
lix = 0
for l in loops:
    print("Loop ", l.id)
    print("  ", l.comps)
for c in comps:
    print(c.name, c.value, c.loops, c.polarity)
rm = makeRarray()
print(rm)
vm = makeVarray()
print(vm)
ires = np.linalg.inv(rm).dot(vm)
print(ires)

get_volt()
get_net_voltages()




