filename = "/home/simao/Sync/article_magnetic/program/kestrel/scripts/tasks.txt"

# Lx_list = []
# Ly_list = []
# seed_list = []
# w_list = []
# nx_list = [] 
# ny_list = []
# B_list = []
# E_list = []
# NR_list = []
# ND_list = []
# N_max_list = []
# N_slice_list = []
# N_current_list = []
# status_list = []

f = open(filename,'r').readlines()
found_next = False
for line in f:
    splitted = line[:-1].split(" ")

    Lx          = splitted[0][3:]
    Ly          = splitted[1][3:]
    seed        = splitted[2][5:]
    w           = splitted[3][2:]
    nx          = splitted[4][3:]
    ny          = splitted[5][3:]
    B           = splitted[6][2:]
    energies    = splitted[7][3:-1]
    NR          = splitted[8][3:]
    ND          = splitted[9][3:]
    N_max       = int(splitted[10][2:])
    N_slice     = int(splitted[11][3:])
    N_current   = int(splitted[12][3:])
    status      = splitted[13][7:]

    # Lx_list.append(Lx)
    # Ly_list.append(Ly)
    # seed_list.append(seed)
    # w_list.append(w)
    # nx_list.append(nx)
    # ny_list.append(ny)
    # B_list.append(B)
    # E_list.append(energies)
    # NR_list.append(NR)
    # ND_list.append(ND)
    # N_max_list.append(N_max)
    # N_slice_list.append(N_slice)
    # N_current_list.append(N_current)
    # status_list.append(status)

    if(status != "R" and status != "C" and not found_next and N_current < N_max):
        print("Lx={0} Ly={1} seed={2} W={3} nx={4} ny={5} B={6} E=[{7}] NR={8} ND={9} N={10} ND={11} NR={12} status={13}".format(Lx, Ly, seed, w, nx, ny, B, energies, NR, ND, N_max, N_slice, N_current, status))
        found_next = True


if not found_next:
    print("done")


# print(Lx_list)
# print(Ly_list)
# print(seed_list)
# print(w_list)
# print(nx_list)
# print(ny_list)
# print(B_list)
# print(E_list)
# print(NR_list)
# print(ND_list)
# print(N_max_list)
# print(N_slice_list)
# print(N_current_list)
# print(status_list)
