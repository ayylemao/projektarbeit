import numpy as np

nat = 8

configs = np.loadtxt("configs")

with open("POSALL.ascii") as f:
    content = f.readlines()


def replaceSymb(iat, configStrings):
    if iat < 4:
        print("replaced C")
        configStrings[iat+2] = configStrings[iat+2].replace("C", "O")
    else:
        print("replaced N")
        configStrings[iat+2] = configStrings[iat+2].replace("N", "F")
    return configStrings


def getConfig(iat, iconfig, nat, content):
    configStrings = []
    iconfig = int(iconfig)
    iat = int(iat)
    for i in range((iconfig-1)*11, (iconfig-1)*11+11):
        configStrings.append(content[i])
    configStrings = replaceSymb(iat, configStrings)
    return configStrings



for i in range(0, 100):
    if i < 10:
        name = "landmark_" + "0" + str(i) + ".ascii"
    else:
        name = "landmark_" + str(i) + ".ascii"
    iat = configs[i][1]
    iconfig = configs[i][2]
    print(iat, iconfig)
    with open(name, 'w') as f:
        strings = getConfig(iat, iconfig, 8, content)
        for i in range(0, 11):
            f.write(strings[i])


# iat = int(configs[0][1])
# iconfig = int(configs[0][2])

# a = getConfig(8, 85, 8, content)

# for i in range(0, 11):
#     print(a[i])



