#!/bin/awk -f

# BEGIN specifies that actions are performed before the fist line of input is read.
# OFS stands for "Output Field Separator" and this code assigns a tab to this value (default is a single space)
BEGIN{
 OFS="\t"
 h["1/1"]=2
 h["0/1"]=1
 h["0/0"]=0
}

# && is a Boolean Operator for Logican AND
# NF stands for "Number of Fields"
# The for loop has 3 parts: the firs is ine initial expression, which assigns the initial value of the counter variable.
# The second part of the loop is a test expression that is evaluated before execution. (If false, the loop is exited).
# The third part increases the counter variable after each pass (++ is an increment operator)
# Patterns contained within forward slashes (such as /^#CHROM/) are specific patterns that AWK should look for to do something in particular.
# Here, the script is told what to expect for the first line of interest in the VCF file and stores it inside s[i] to be printed later as output.

/^#CHROM/ && !header {
 for(i = 1; i <= NF; i++)
  s[i]=$i
 header=1
}

# || is a Boolean Operator for Logical OR
# These lines cause awk to discard the current input record and proceed to the next
# ~ is a relational expression that means "matches regular expression"
/INDEL/ || $4 == "N" || $5 ~ /,/ || /^#/ { 
 next
}


# Another for loop.
# split() splits a string into an array (called "a" here) which is delimited by ":"
# Output prints contents of fields (columns) 1,2,4,5, and 6 plus header information from s[i] and
{
 for(i = 10; i <= NF; i++){
  split($i, a, ":")
  if(a[3] > 0)
   print $1, $2, $4, $5, $6, s[i], h[a[1]], a[5], a[3], a[4]
 }
}
