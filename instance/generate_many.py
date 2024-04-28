from generate import gen_one

def main():
	for (n1, n2) in [
				(4, 4),
				(5, 5),
				(6, 6),
				(6, 7),
				(7, 7),
				(7, 8),
				(8, 8),
			]:
		for d in [20, 30, 40, 50, 60, 70, 80]:
			gen_one(n1, n2, d)


if __name__ == "__main__":
    main()
