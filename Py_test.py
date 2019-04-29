###
### Adding Eight Numbers
### Comp. 11 HW1 
### Author: Will Kendall 1/28/18
###

#
#User input
#
numInput = input("Eww! I'm drowning I need exactly eight floats to survive: ")

#
#Splits string by " "
#
num_data = numInput.split()

#
#Nested functions: The map functions applys float to all strings in 
#num_data, the list function lists those floats
#
num_datafin = list(map(float, num_data))

#
#Finds length of floated number list
#
num_datalen = len(num_datafin)

#
#Conditional: Does not print sum unless user inputs 8 numbers. Elifs 
#provide alternates.
#

if num_datalen < 8:
  print('*glug glug glug...*')
elif num_datalen > 8:
  print("Too many floats! I'm floating awayyyyyyyy")
else:
    print('You saved me brave hero. Your reward is the sum of your floats: ')
    numsum = sum(num_datafin)
    print(numsum, 'Coins')


