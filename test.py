
string = "234–133–96–64–83–184–179–168–84–248–125–147–155–167–215–288–160–252–151–322– 128–210–144–307–185–151–69–167–131–197–157–127–277–252–126–340–241–140–70–123– 85–168–224–126–100–249–82–165–111–133–252–287–184–288–292–259–92–82–221–308–313– 182–88–374–129–169–231–239–343–396–194–340–337–336–80–211–275–386–382–192–116– 387–122–189–267–235–464–459–216–454"

count = 1

for char in string:
    if char == "–":
        count += 1


print(count)