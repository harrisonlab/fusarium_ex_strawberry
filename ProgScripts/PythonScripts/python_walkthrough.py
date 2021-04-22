# Game Over - Version 2
# Demonstrates the use of quotes in strings

# print - echos what you type verbatum
# need to use brackets and "" to verify
print("Program 'Game Over' 2.0")
print("Same", "message", "as before")
print("Just",
        "a bit",
        "bigger")

print("Here", end=" ")
print("it is...")
# Triple-quoted strings span multiple lines. They print exactly the way you type them.
print(
        """
        _____ ___ ___ ___ _____
        / ___| / | / |/ | | ___|
        | | / /| | / /| /| | | |__
        | | _ / ___ | / / |__/ | | | __|
        | |_| | / / | | / / | | | |___
        \_____/ /_/ |_| /_/ |_| |_____|
        _____ _ _ _____ _____
        / _ \ | | / / | ___| | _ \
        | | | | | | / / | |__ | |_| |
        | | | | | | / / | __| | _ /
        | |_| | | |/ / | |___ | | \ \
        \_____/ |___/ |_____| |_| \_\

        """
        )

# Used the tab escape sequence, \t, three times in a row. So, when it prints the string, it prints three tabs and then fancy credits
print("\t\t\tFancy Credits")

print("\t\t\t \\ \\ \\ \\ \\ \\ \\")
print("\t\t\t\tby")
print("\t\t\tMichael Dawson")
print("\t\t\t \\ \\ \\ \\ \\ \\ \\")

print("\nSpecial thanks goes out to:")
print("My hair stylist, Henry \'The Great,\' who never says \"can\'t.\"")

print("\a")

input("\n\nPress the enter key to exit.")
# input - a function allows you to input a step control by fulfilling a condition

# print ice cream flavour

print("Vanilla")

print("Injustice anywhere is a threat to justice everywhere.")
print("MLK Jr.")

input("\n\nPress the enter key to exit.")

# Word Problems
# Demonstrates numbers and math

# set a variable - a = 16
# putting an '=' after the math function and before the value you are multiplying changes the value of your variable
#  a *=4 - this means a would then become 16 * 4 which is 64

# abs - absolute distance from zero
# abs(5) is 5
# abs(-5) is also 5

print("If a 2000 pound pregnant hippo gives birth to a 100 pound calf,")
print("but then eats 50 pounds of food, how much does she weigh?")
input("Press the enter key to find out.")
print("2000 - 100 + 50 =", 2000 - 100 + 50)
# + - are same symbols

print("\nIf an adventurer returns from a successful quest and buys each of")
print("6 companions 3 bottles of ale, how many bottles are purchased?")
input("Press the enter key to find out.")
print("6 * 3 =", 6 * 3)
# * is multiply

print("\nIf a restaurant check comes to 19 dollars with tip, and you and")
print("your friends split it evenly 4 ways, how much do you each throw in?")
input("Press the enter key to find out.")
print("19 / 4 =", 19 / 4)
# forward slash is divide

print("\nIf a group of 4 pirates finds a chest full of 107 gold coins, and")
print("they divide the booty evenly, how many whole coins does each get?")
input("Press the enter key to find out.")
print("107 // 4 =", 107 // 4)
# double forward slash is whole number division  (integer)

print("\nIf that same group of 4 pirates evenly divides the chest full")
print("of 107 gold coins, how many coins are left over?")
input("Press the enter key to find out.")
print("107 % 4 =", 107 % 4)
input("\n\nPress the enter key to exit.")
# % gives remainder (modulus)

print("\nIf a ball has a diameter of 20cm,")
print("what is the volume of the ball?")
print("\nThe volume of a sphere is 4/3 * pi * 10e3")
print("(4 /3 ) * 3.14 * 10 **3=," (4 /3 ) * 3.14 * 10**3)
input("\n\nPress the enter key to exit.")
# 10**3 - the to ** indicate to the power of

name = "Akinkunmi"
# x = abc - setting a variable
# sets value for variable - don't need $ like in linux

print(name)
print("Hi,", name)
input("\n\nPress the enter key to exit.")

# Makes variable an input decided by user
name1 = input("Hi. What's your name? ")

print(name1)

print("Hi,", name1)

input("\n\nPress the enter key to exit.")

# Quotation Manipulation
# Demonstrates string methods

# quote from IBM Chairman, Thomas Watson, in 1943
quote = "I think there is a pokemon under my bed."
print("Original quote:")
print(quote)

#quote.upper invokes argument inside function
print("\nIn uppercase:")
print(quote.upper())

print("\nIn lowercase:")
print(quote.lower())

print("\nAs a title:")
print(quote.title())

# prints a new string with replacement
print("\nWith a minor replacement:")
print(quote.replace("under", "on"))

print("\nOriginal quote is still:")
print(quote)

print("\nWith swapped cases:")
print(quote.swapcase())

#Returns a new string where the first letter is capitalized and the rest are lowercase.
print("\nWith capitalS:")
print(quote.capitalize())

input("\n\nPress the enter key to exit.")

# Trust Fund Buddy - Good
# Demonstrates type conversion

print(
"""
                        Trust Fund Buddy

Totals your monthly spending so that your trust fund doesn't run out
(and you're forced to get a real job).

Please enter the requested, monthly costs. Since you're rich, ignore
pennies and use only dollar amounts.

"""
)
# Making each input an integer allows total to have numerical value
# variable = int(input("abc:"))

car = input("Lamborghini Tune-Ups: ")
car = int(car)

rent = int(input("Manhattan Apartment: "))
jet = int(input("Private Jet Rental: "))
gifts = int(input("Gifts: "))
food = int(input("Dining Out: "))
staff = int(input("Staff (butlers, chef, driver, assistant): "))
guru = int(input("Personal Guru and Coach: ") )
games = int(input("Computer Games: "))

total = car + rent + jet + gifts + food + staff + guru + games

print("\nGrand Total:", total)
input("\n\nPress the enter key to exit.")
